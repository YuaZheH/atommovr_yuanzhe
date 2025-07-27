"""A customized Gymnasium environment for atom rearrangement. 
Designed to be a parent class for all rearrangement environments.

Implements the key elements of a Gymnsium environment:
1. __init__() (the structure for the action space and observation space)
2. step() (the function which applies a DRL agent's chosen action to the current observation)
3. reset()
4. render() (placeholder)
    
NB: this has many placeholder functions (e.g. get_reward() is not defined), so it is NOT suitable
to be run on its own. Instead, you should build a child class from it.

MO: we use the `gymnasium module`` to create the "game environment" 
that specifies the rearrangement state and possible actions. """

import numpy as np
import gymnasium as gym
from typing import Any

# to set rearrangement rules and make visualizations
from atommover.utils.core import PhysicalParams
from atommover.utils.move_utils import move_atoms, get_move_list_from_AOD_actions

class ParentRearrangementEnv(gym.Env):
    """ This is a parent class for all rearrangement environments. 
        It contains mostly placeholder functions and should not be used for training or testing. 
        
        There are four key functions in every Gymnasium environment:
        1. `__init__()`
            
            This defines the key elements of the environment, namely the action space, the observation space, and the current state.
        2. `step()`
            
            This takes as input the agent's chosen action and applies it to the current state, returning a new state.
        3. `render()`
            
            We ignore this for the time being, but it's used for generating a visualization of the current state. For us, this is generally handled outside of the environment (see the function `visualize_testing_episode` in `source/drl/train_funcs.py`).
        4. `reset()`
            
            This is an operation to restart the environment after the end of a training episode.  """
    def __init__(self, target_config, params: PhysicalParams = PhysicalParams()):
        super(ParentRearrangementEnv, self).__init__()
        # Space of actions we can take
        # these correspond to the horizontal and vertical AOD frequencies
        self.array_len_x = len(target_config[0])
        self.array_len_y = len(target_config)
        AOD_action_shapes = 4*np.ones(self.array_len_x+self.array_len_y, dtype = int)
        # action space is conceptualized as array_len_x + array_len_y discrete action spaces, representing AOD frequencies:
        # each action space is discrete-4: NOOP/OFF[0], TURNON[1], RAMPUP[2], RAMPDOWN[3]
        self.action_space = gym.spaces.MultiDiscrete(AOD_action_shapes, start = np.zeros(self.array_len_x+self.array_len_y))

        # Observation space
        self.observation_space = gym.spaces.Box(low=0, high=1, shape = (2,self.array_len_y, self.array_len_x), dtype = np.uint8)
        
        # Set parameters
        self.params = params

        # Generate initial and target states randomly, with the desired number of atoms
        self.state = np.zeros((2, self.array_len_y, self.array_len_x), dtype = np.uint8)
        self.target = self._get_init_state()
        init_config = self._get_init_state()
        self.state[0,:,:] = init_config
        self.state[1,:,:] = self.target

        # Set time
        self.time = 0

    def _get_init_state(self):
        pass

    def get_reward(self,observation, delta_time):
        pass

    def step(self, action):
        # breaking action into AOD commands
        horizontal_AOD_cmds = action[0:self.array_len_x]
        vertical_AOD_cmds = action[self.array_len_x:(self.array_len_y+self.array_len_x)]

        # getting move list from AOD commmands
        move_list = get_move_list_from_AOD_actions(horizontal_AOD_cmds, vertical_AOD_cmds)
        if np.sum(self.state[1,:,:]) >= len(move_list) > 0: # NKH 01/03/25: made this so it wouldn't make unnecessary moves. OG: # if len(move_list) > 0:
            print(len(move_list))
            # running moves in errorless environment
            observation, _ = move_atoms(self.state[0,:,:].reshape(self.array_len_y, self.array_len_x), 
                                        move_list, 
                                        look_for_flag = True)
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
            surplus_atoms = int(np.sum(self.state[0,self.padding_y0:(self.array_len_y-self.padding_y1),self.padding_x0:(self.array_len_x-self.padding_x1)])) - int(np.sum(self.state[1,self.padding_y0:(self.array_len_y-self.padding_y1),self.padding_x0:(self.array_len_x-self.padding_x1)]))
        except NameError:
            print('`self.padding_y0` or `self.padding_x0` does not exist')
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
    
    def render(self):
        pass

    def reset(self, seed = None, options: 'dict[str, Any]' = None):
        super().reset(seed=seed, options=options)
        self.target = self._get_init_state()
        init_is_different = False
        while not init_is_different:
            init_config = self._get_init_state()
            if not np.array_equal(init_config, self.target):
                init_is_different = True

        self.state[0,:,:] = init_config
        self.state[1,:,:] = self.target
        self.time = 0
        return self.state, {}