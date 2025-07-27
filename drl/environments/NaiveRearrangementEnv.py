# This is a naive example of how to build a complete environment that you can use for training,
# using ParentRearrangementEnv as a parent class.

import numpy as np
import gymnasium as gym

# to compute loss function
from ot.coot import co_optimal_transport2 as coot2

# to set rearrangement rules and make visualizations
from atommover.utils.core import random_loading
from drl.environments.ParentRearrangementEnv import ParentRearrangementEnv

class NaiveRearrangementEnv(ParentRearrangementEnv):
    """ This is a naive example of how to build a complete environment that you can use for training. 
        It takes `ParentRearrangementEnv` as a parent class."""

    def get_init_state(self):
        init_config = random_loading([self.array_len_x, self.array_len_y], self.params.loading_prob)
        return init_config

    def get_reward(self,observation):
        init_distance = self.array_len_x*self.array_len_y*coot2(self.state,self.target)
        final_distance = self.array_len_x*self.array_len_y*coot2(observation,self.target)
        reward = init_distance-final_distance - self.time*100
        if np.array_equal(observation, self.target): #TODO: fix this (penalizes having extra atoms)
            reward += 50
        return reward