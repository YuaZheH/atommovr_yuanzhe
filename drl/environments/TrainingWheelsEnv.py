# An environment for easy training. 

import numpy as np
import gymnasium as gym
from typing import Any

# to set rearrangement rules and make visualizations
from atommover.utils.core import PhysicalParams, random_loading
from atommover.utils.move_utils import move_atoms, get_move_list_from_AOD_actions
from drl.environments.ParentRearrangementEnv import ParentRearrangementEnv

class TrainingWheelsEnv(ParentRearrangementEnv):
    """ 
    Builds off of `ParentRearrangementEnv` to create an environment with tunable difficulty for training the RL agent.
    
    The difficulty of the environment can be increased by increasing `n_atoms` and decreasing `padding`.
    
    At the start of each episode, prepares a random initial configuration and a random target configuration for the agent to prepare.
    
    ### Args: 
    `size` (list of ints): 
        the number of rows and columns in the array. Must be 2D.
    
    `params` (`PhysicalParams()` object): 
        specifies the atom, tweezer, and array properties.
    
    `padding` (int): 
        the number of rows or columns to ignore, counted from each edge of the array.
    
    `n_atoms` (int): 
        the number of atoms in the target array.
    
    `max_actions_per_episode` (int): 
        the maximum number of actions the RL agent can take to prepare a given configuration. 
        Used to force the agent to find fast solutions. """
    def __init__(self, size: list = [11,13], params = PhysicalParams(), padding: list = [[5, 5],[3, 3]], n_atoms = 1, max_actions_per_episode = 5) -> None:
        # Esssential Parameters
        self.array_len_x = size[1] # the number of columns in the array
        self.array_len_y = size[0] # the number of rows in the array
        self.padding_x0, self.padding_x1 = padding[1] # the number of columns to disable actions for (from each edge)
        self.padding_y0, self.padding_y1 = padding[0] # the number of rows to disable actions for (from each edge)
        self.n_atoms = n_atoms # the number of atoms in the array
        self.max_actions_per_episode = max_actions_per_episode # the maximum number of actions the agent can take per training episode (used to force the agent to learn fast)
        if self.array_len_x <= (self.padding_x0+self.padding_x1) or self.array_len_y <= (self.padding_y0+self.padding_y1):
            raise Exception(f'Padding of {padding} is too large for array of size ({self.array_len_x},{self.array_len_y}).')
        
        # Action space
            # these correspond to the horizontal and vertical AOD frequencies
        AOD_action_shapes = 4*np.ones(self.array_len_x+self.array_len_y, dtype = int)

            # action space is conceptualized as array_len_x + array_len_y discrete action spaces, representing AOD frequencies:
            # each action space is discrete-4: NOOP/OFF[0], TURNON[1], RAMPUP[2], RAMPDOWN[3]
        self.action_space = gym.spaces.MultiDiscrete(AOD_action_shapes, start = np.zeros(self.array_len_x+self.array_len_y))

        # Observation space
        self.observation_space = gym.spaces.Box(low=0, high=1, shape = (2, self.array_len_y, self.array_len_x), dtype = np.uint8)
        
        # Set experimental parameters
        self.params = params

        # Generate initial and target states randomly, with the desired number of atoms
        self.state = np.zeros((2, self.array_len_y, self.array_len_x), dtype = np.uint8)
        target = self._get_init_state()
            # we make sure that the configuration we prepare is different from the target configuration
        init_is_different = False
        while not init_is_different:
            init_config = self._get_init_state()
            if not np.array_equal(init_config, target):
                init_is_different = True

        self.state[0,:,:] = init_config
        self.state[1,:,:] = target

        # Set time
        self.time = 0

    def _get_init_state(self):
        """ This function prepares a random array of atoms with exactly `self.n_atoms` atoms. """
        sufficient = False
        effective_len_x = self.array_len_x-self.padding_x0-self.padding_x1
        effective_len_y = self.array_len_y-self.padding_y0-self.padding_y1
        while sufficient == False:
            init_config = np.zeros([1,self.array_len_y, self.array_len_x])
            loading_prob = self.n_atoms/(effective_len_x*effective_len_y)
            middle = random_loading([effective_len_y, effective_len_x], loading_prob)
            init_config[0,self.padding_y0:(self.array_len_y-self.padding_y1),self.padding_x0:(self.array_len_x-self.padding_x1)] = middle

            if np.sum(init_config) == self.n_atoms:
                sufficient = True
        return init_config

    def get_reward(self,observation, delta_time) -> float:
        # init_distance = self.array_len_x*self.array_len_y*coot2(self.state.reshape(self.array_len_x, self.array_len_y),
        #                                                         self.target.reshape(self.array_len_x, self.array_len_y))
        # final_distance = self.array_len_x*self.array_len_y*coot2(observation.reshape(self.array_len_x, self.array_len_y),
        #                                                          self.target.reshape(self.array_len_x, self.array_len_y))
        # reward = init_distance-final_distance - delta_time*100

        if np.array_equal(observation.reshape(self.array_len_y, self.array_len_x), self.state[1,:,:].reshape(self.array_len_y, self.array_len_x)):
            reward = 50.0
        elif np.sum(observation[self.padding_y0:self.array_len_y-self.padding_y1, self.padding_x0:self.array_len_x-self.padding_x1]) < np.sum(self.state[1][self.padding_y0:self.array_len_y-self.padding_y1, self.padding_x0:self.array_len_x-self.padding_x1]):
            reward = -50.0
        else:
            reward = -100*delta_time
        return reward
    
    def get_action_mask(self) -> 'tuple[np.ndarray]':
        """ This code generates action masks for arrays with padding 
            (i.e. arrays where we only want to operate on some subarray).
            
            This type of action mask is intended for training purposes."""
        
        # defining base masks
        NOOP_mask = np.array([1,0,0,0], dtype=np.int8)
        ALLOP_mask = np.array([1,1,1,1], dtype=np.int8)
        mask_list = []

        ### Setting masks for horizontal/col operations ##
            # for cols where there are no atoms, disable atom moves
        for _ in range(self.padding_x0):
            mask_list.append(NOOP_mask)
        
            # leftmost col
        if self.array_len_x-self.padding_x0-self.padding_x1-2 >= 0:
            mask_list.append(np.array([1,1,1,0], dtype=np.int8))
        else:
            mask_list.append(np.array([0,1,0,0], dtype=np.int8))

            # middle cols
        for _ in range(self.array_len_x-self.padding_x0-self.padding_x1-2):
            mask_list.append(ALLOP_mask)

            # rightmost col
        mask_list.append(np.array([1,1,0,1], dtype=np.int8))
        
            # for cols where there are no atoms, disable atom moves
        for _ in range(self.padding_x1):
            mask_list.append(NOOP_mask)
        
        ### Setting masks for vertical/row operations ##
            # for rows where there are no atoms, disable atom moves
        for _ in range(self.padding_y0):
            mask_list.append(NOOP_mask)
        
            # topmost row
        if self.array_len_y-self.padding_y0-self.padding_y1-2 >= 0:
            mask_list.append(np.array([1,1,1,0], dtype=np.int8))
        else:
            mask_list.append(np.array([0,1,0,0], dtype=np.int8))

            # middle rows
        for _ in range(self.array_len_y-self.padding_y0-self.padding_y1-2):
            mask_list.append(ALLOP_mask)

            # bottommost row #TODO check this
        if self.array_len_y-self.padding_y0-self.padding_y1-2 >= 0:
            mask_list.append(np.array([1,1,0,1], dtype=np.int8))
        
            # for rows where there are no atoms, disable atom moves
        for _ in range(self.padding_y1):
            mask_list.append(NOOP_mask)

        return tuple(mask_list)
    

class BabyEnv(TrainingWheelsEnv):
    """ 
    Builds off of `ParentRearrangementEnv` to create an environment with tunable difficulty for training the RL agent.
    
    The difficulty of the environment can be increased by increasing `n_atoms` and decreasing `padding`.
    
    At the start of each episode, prepares a random initial configuration and a random target configuration for the agent to prepare.
    
    ### Args: 
    `size` (list of ints): 
        the number of rows and columns in the array. Must be 2D.
    
    `params` (`PhysicalParams()` object): 
        specifies the atom, tweezer, and array properties.
    
    `padding` (int): 
        the number of rows or columns to ignore, counted from each edge of the array.
    
    `n_atoms` (int): 
        the number of atoms in the target array.
    
    `max_actions_per_episode` (int): 
        the maximum number of actions the RL agent can take to prepare a given configuration. 
        Used to force the agent to find fast solutions. """
    def __init__(self, size: list, params = PhysicalParams(), padding = list, n_atoms = 1, max_actions_per_episode = 5) -> None:
        # Esssential Parameters
        self.array_len_x = size[1] # the number of columns in the array
        self.array_len_y = size[0] # the number of rows in the array
        self.padding_x = padding[1] # the number of columns to disable actions for (from each edge)
        self.padding_y = padding[0] # the number of rows to disable actions for (from each edge)
        self.n_atoms = n_atoms # the number of atoms in the array
        self.max_actions_per_episode = max_actions_per_episode # the maximum number of actions the agent can take per training episode (used to force the agent to learn fast)
        if self.array_len_x <= self.padding_x*2 or self.array_len_y <= self.padding_y*2:
            raise Exception(f'Padding of {padding} is too large for array of size ({self.array_len_x},{self.array_len_y}).')
        
        # Action space
            # these correspond to the horizontal and vertical AOD frequencies
        AOD_action_shapes = 4*np.ones(self.array_len_x+self.array_len_y, dtype = int)

            # action space is conceptualized as array_len_x + array_len_y discrete action spaces, representing AOD frequencies:
            # each action space is discrete-4: NOOP/OFF[0], TURNON[1], RAMPUP[2], RAMPDOWN[3]
        self.action_space = gym.spaces.Discrete(3)

        # Observation space
        self.observation_space = gym.spaces.Box(low=0, high=1, shape = (2, self.array_len_y, self.array_len_x), dtype = np.uint8)
        
        # Set experimental parameters
        self.params = params

        # Generate initial and target states randomly, with the desired number of atoms
        self.state = np.zeros((2, self.array_len_y, self.array_len_x), dtype = np.uint8)
        target = self._get_init_state()
            # we make sure that the configuration we prepare is different from the target configuration
        init_is_different = False
        while not init_is_different:
            init_config = self._get_init_state()
            if not np.array_equal(init_config, target):
                init_is_different = True
        self.state[0,:,:] = init_config
        self.state[1,:,:] = target

        # Set time
        self.time = 0

    def get_action_mask(self) -> 'tuple[np.ndarray]':
        """ This code generates action masks for arrays with padding 
            (i.e. arrays where we only want to operate on some subarray).
            
            This type of action mask is intended for training purposes."""
        
        # defining base masks
        ALLOP_mask = np.array([1,1,1], dtype=np.int8)
        mask_list = [ALLOP_mask]

        return tuple(mask_list)
    
    def step(self, action):
        # breaking action into AOD commands
        vertical_AOD_cmds = [0,0,0,0,0,1,0,0,0,0,0]
        if action == 0:
            # do nothing
            horizontal_AOD_cmds = [0,0,0,0,0,0,0,0,0,0]
        elif action == 1:
            # move atom to right
            horizontal_AOD_cmds = [0,0,0,0,2,0,0,0,0,0]
        elif action == 2:
            # move atom to left
            horizontal_AOD_cmds = [0,0,0,0,0,3,0,0,0,0]
        # getting move list from AOD commmands
        move_list = get_move_list_from_AOD_actions(horizontal_AOD_cmds, vertical_AOD_cmds)
        if len(move_list) > 0:
            # running moves in errorless environment
            observation, _ = move_atoms(self.state[0,:,:], 
                                        move_list, 
                                        look_for_flag = False, 
                                        putdown_fail_rate=self.params.putdown_fail_rate, 
                                        pickup_fail_rate=self.params.pickup_fail_rate)
        else:
            observation = self.state[0,:,:]

        # adding time step
        delta_time = self.params.spacing/self.params.AOD_speed
        self.time += delta_time

        # cost function for calculating reward
        reward = self.get_reward(observation, delta_time)

        # update current state (NB: must be AFTER calculating reward)
        self.state[0,:,:] = observation.reshape(1,self.array_len_y, self.array_len_x)
        
        # determine whether the episode terminates
        terminated = np.array_equal(self.state[0,:,:], self.state[1,:,:])
        try:
            surplus_atoms = int(np.sum(self.state[0,self.padding_y:(self.array_len_y-self.padding_y),self.padding_x:(self.array_len_x-self.padding_x)])) - int(np.sum(self.state[1,self.padding_y:(self.array_len_y-self.padding_y),self.padding_x:(self.array_len_x-self.padding_x)]))
        except NameError:
            print('`self.padding_y` or `self.padding_x` does not exist')
            surplus_atoms = int(np.sum(self.state[0,:,:])) - int(np.sum(self.state[1,:,:]))
        if surplus_atoms < 0:
            terminated = True
            reward += surplus_atoms

        try:
            if self.max_actions_per_episode != 0 and self.time/delta_time >= self.max_actions_per_episode:
                terminated = True
                reward -= 50
        except AttributeError:
            pass

        # placeholders for `info` and `truncated` (gym.Env parameters that we don't use)
        info = {}
        truncated = False

        return self.state, reward, terminated, truncated, info

