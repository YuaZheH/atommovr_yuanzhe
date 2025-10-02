# Single-species algorithms. 

# FOR CONTRIBUTORS:
# - Please write your algorithm in a separate .py file
# - Once you have done that, please make an algorithm class with the following three functions (see the `Algorithm` class for more details)
#   1. __repr__(self) - this should return the name of your algorithm, to be used in plots.
#   2. get_moves(self) - given an AtomArray object, returns a list of lists of Move() objects.
#   3. (optional) __init__() - if your algorithm needs to use arguments that cannot be specified in AtomArray
import numpy as np

from atommover.utils.AtomArray import AtomArray
from atommover.algorithms.Algorithm_class import Algorithm
from atommover.algorithms.source.balance_compact import balance_and_compact
from atommover.algorithms.source.bc_new import bcv2
from atommover.algorithms.source.generalized_balance import generalized_balance
# from atommover.algorithms.source.Hungarian_works import parallel_Hungarian_algorithm_works, parallel_LBAP_algorithm_works, Hungarian_algorithm_works
from atommover.algorithms.source.Hungarian_works import parallel_Hungarian_algorithm_works, Hungarian_algorithm_works


##########################
# Bernien Lab algorithms #
##########################

# Parallel Hungarian
class ParallelHungarian(Algorithm):
    """ A variant on the Hungarian matching algorithm that parallelizes the moves
        instead of executing them sequentially (one by one).
        
        Supported configurations: all. """
    def __repr__(self):
        return 'Parallel Hungarian'
    
    def get_moves(self, atom_array: AtomArray, do_ejection: bool = False, final_size: list = [], round_lim: int = 0):
        if atom_array.n_species != 1:
            raise ValueError(f"Single-species algorithm cannot process atom array with {atom_array.n_species} species.")
        if round_lim == 0:
            round_lim = int(np.sum(atom_array.target))
        return parallel_Hungarian_algorithm_works(atom_array.matrix, atom_array.target, do_ejection, final_size, round_lim)

# class ParallelLBAP(Algorithm):
#     """ Solves the linear bottleneck assignment problem and parallelizes the moves.
#         Code taken from ParallelHungarian.
        
#         Supported configurations: all. """
#     def __repr__(self):
#         return 'Parallel LBAP'
    
#     def get_moves(self, atom_array: AtomArray, do_ejection: bool = False, final_size: list = [], round_lim: int = 0):
#         if atom_array.n_species != 1:
#             raise ValueError(f"Single-species algorithm cannot process atom array with {atom_array.n_species} species.")
#         if round_lim == 0:
#             round_lim = int(np.sum(atom_array.target))
#         return parallel_LBAP_algorithm_works(atom_array.matrix, atom_array.target, do_ejection, round_lim)

# Generalized Balance
class GeneralizedBalance(Algorithm):
    """Implements the generalized balance algorithm, which alternatively operates 
       row balance and column balance algorithms, as originally described by Bo-Yu
       and Nikhil in the Bernien lab meeting GM 268.
       
       Supported configurations: all. """

    def __repr__(self):
        return 'Generalized Balance'
    
    def get_moves(self, atom_array: AtomArray, do_ejection: bool = False):
        if atom_array.n_species != 1:
            raise ValueError(f"Single-species algorithm cannot process atom array with {atom_array.n_species} species.")
        return generalized_balance(atom_array.matrix[:,:,0], atom_array.target, do_ejection)

###########################################
# Existing algorithms from the literature #
###########################################

# Hungarian
class Hungarian(Algorithm):
    """ Implements the Hungarian matching algorithm, which generates a cost 
        matrix mapping available atoms to the target spots, and solves the
        linear assignment problem to find an efficient set of moves.
        
        Supported configurations: all. """

    def __repr__(self):
        return 'Hungarian'

    def get_moves(self, atom_array: AtomArray, do_ejection: bool = False):
        if atom_array.n_species != 1:
            raise ValueError(f"Single-species algorithm cannot process atom array with {atom_array.n_species} species.")
        return Hungarian_algorithm_works(atom_array.matrix[:,:,0], atom_array.target, do_ejection)

# Balance and Compact
class BCv2(Algorithm):
    """Implements the Balance and Compact algorithm, as originally described
       in [PRA 70, 040302(R) (2004)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.70.040302)
       
       Supported configurations: `Configurations.MIDDLE_FILL`"""

    def __repr__(self):
        return 'Balance & Compact'

    def get_moves(self, atom_array: AtomArray, do_ejection: bool = False):
        if atom_array.n_species != 1:
            raise ValueError(f"Single-species algorithm cannot process atom array with {atom_array.n_species} species.")
        return bcv2(atom_array, do_ejection)
    


# Balance and Compact
class BalanceAndCompact(Algorithm):
    """ NOTE: we recommend that you use the (faster) BCv2 algorithm. 
        This is an older version that we used to make Fig. 2 in the paper.

        A slow implementation of the Balance and Compact algorithm, as originally described
        in [PRA 70, 040302(R) (2004)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.70.040302)
       
        Supported configurations: `Configurations.MIDDLE_FILL`"""

    def __repr__(self):
        return 'Balance & Compact (slow)'

    def get_moves(self, atom_array: AtomArray, do_ejection: bool = False):
        if atom_array.n_species != 1:
            raise ValueError(f"Single-species algorithm cannot process atom array with {atom_array.n_species} species.")
        return balance_and_compact(atom_array.matrix[:,:,0], atom_array.target, do_ejection)
