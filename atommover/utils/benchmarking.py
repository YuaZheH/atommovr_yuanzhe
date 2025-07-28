# Object for running benchmarking rounds and saving data


import sys
import copy
import random
import math
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from atommover.utils.errormodels import ZeroNoise
from atommover.utils.core import generate_random_target_configs, generate_random_init_configs, PhysicalParams, Configurations, CONFIGURATION_PLOT_LABELS
from atommover.utils.move_utils import move_atoms
from atommover.utils.AtomArray import AtomArray
from atommover.algorithms.Algorithm_class import Algorithm, get_effective_target_grid

# ## Sanity check ##
# def evaluate_moves_old(init_config: np.ndarray, move_list: list, params: PhysicalParams):
#     # making reference time
#     t_total = 0
#     N_parallel_moves = 0
#     N_non_parallel_moves = 0
#     matrix = np.array(copy.deepcopy(init_config))

#     # iterating through moves and updating matrix
#     for move_ind, move_set in enumerate(move_list):

#         # performing the move
#         matrix, [failed_moves, flags] = move_atoms(matrix, move_set, look_for_flag = True, pickup_fail_rate = params.pickup_fail_rate, putdown_fail_rate = params.putdown_fail_rate)
#         N_parallel_moves += 1
#         N_non_parallel_moves += len(move_set)

#         # calculating the distance of each move in the set
#         distances = []
#         for move_set_ind, move in enumerate(move_set):
#             distances.append(move.distance)

#         # calculating the time to complete the move set in parallel
#         t_move = max(distances)*params.spacing/params.AOD_speed
#         t_total += t_move

#         # simulating atom loss
#         matrix, loss_flag = atom_loss(matrix, t_move, params.lifetime)

#     return matrix, float(t_total), [N_parallel_moves, N_non_parallel_moves]

def evaluate_moves(array: AtomArray ,move_list: list):
    # making reference time
    t_total = 0
    N_parallel_moves = 0
    N_non_parallel_moves = 0

    # iterating through moves and updating matrix
    for move_ind, move_set in enumerate(move_list):

        # performing the move
        [failed_moves, flags], move_time = array.move_atoms(move_set)
        N_parallel_moves += 1
        N_non_parallel_moves += len(move_set)

        # calculating the time to complete the move set in parallel
        t_total += move_time

    return array, float(t_total), [N_parallel_moves, N_non_parallel_moves]


class BenchmarkingFigure():
    """ 
        NB: this No longer supported.
    Class that specifies plot parameters and figure types to be used in conjunction with the `Benchmarking` class.
    
    NB: this class just specifies what you want to plot, to actually plot you have to pass it to an instance of the 
    `Benchmarking` class and call the `plot_results()` function.

    ## Parameters
    - `y_axis_variables` (list):
        the observables to plot. Must be in ['Success rate', 'Filling fraction', 'Time', 'Wrong places #', 'Total atoms']
    - `figure_type` (str):
        The kind of figure you want to make. Options are histogram ('hist'), a plot comparing different algorithms ('scale'), or a plot comparing different target configurations for the same algorithm ('pattern').
    """
    def __init__(self, variables: list = ['Success rate'], figure_type: str = 'scale'):
        for variable in variables:
            if variable not in ['Success rate', 'Filling fraction', 'Time', 'Wrong places #', 'Total atoms']:
                raise KeyError(f"Variable '{variable}' is not recognized. The only allowed variables are the following: ['Success rate', 'Filling fraction', 'time', 'Wrong places #', 'Total atoms'].")
        self.y_axis_variables = variables
        self.figure_type = figure_type

    def generate_scaling_figure(self, x_axis, benchmarking_results, title, x_label, save, savename = 'Algorithm_scaling'):

        # Iterate over the y-axis variables
        fig, ax = plt.subplots(len(self.y_axis_variables), 1, figsize = (5, 5*len(self.y_axis_variables)))
        for varind, y_var in enumerate(self.y_axis_variables):
            n_datapoints_added = 0
            y_axis = []

            # Iterate over the benchmarking results of each algorithm
            for algo_results in benchmarking_results:
                # If the y-axis variable is a list (e.g. filling fraction), take its average
                if type(algo_results[y_var]) is list:
                    algo_results[y_var] = np.mean(algo_results[y_var])

                if math.isnan(algo_results[y_var]):
                    raise Exception("Data to plot contains nan, indicating that something went wrong in your benchmarking. Please examine data and try again.")
                y_axis.append(algo_results[y_var])
                
                n_datapoints_added += 1

                # If all the results of the algorithm are collected, plot the results
                if n_datapoints_added % len(x_axis) == 0:
                    try:
                        ax[varind].scatter(x_axis, y_axis, marker = 'o', label = algo_results["algorithm"].__class__.__name__)
                    except TypeError:
                        ax.scatter(x_axis, y_axis, marker = 'o', label = algo_results["algorithm"].__class__.__name__)
                    y_axis = []

            try:
                ax[varind].set_xlabel(x_label)
                ax[varind].set_ylabel(y_var.capitalize())
                ax[varind].set_title(title)
                ax[varind].legend(loc='best')
            except TypeError:
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_var.capitalize())
                ax.set_title(title)
                ax.legend(loc='best')

        if save:
            plt.savefig(f'./figs/'+savename)

    def generate_histogram_figure(self, benchmarking_results, title, x_label, save = False, savename = 'Histogram'):
        hist_data = []
        algos_name = []
        fig, ax = plt.subplots(len(self.y_axis_variables), 1, figsize = (5, 5*len(self.y_axis_variables)))
        for varind, y_var in enumerate(self.y_axis_variables):

            for algo_results in benchmarking_results:
                hist_data.append(algo_results[y_var])
                algos_name.append(str(algo_results["algorithm"]))

            try:
                ax[varind].set_xlabel(y_var.capitalize())
                ax[varind].set_ylabel('Frequency')
                ax[varind].set_title(f"{y_var.capitalize()} histogram")
                ax[varind].hist(hist_data, bins = 10, label = algos_name)
                ax[varind].legend()
            except TypeError:
                ax.set_xlabel(y_var.capitalize())
                ax.set_ylabel('Frequency')
                ax.set_title(f"{y_var.capitalize()} histogram")
                ax.hist(hist_data, bins = 10, label = algos_name)
                ax.legend()
        
        if save:
            plt.savefig(f'./figs/{savename}')

    
    def generate_pattern_figure(self, x_axis, benchmarking_results, title, x_label, save = False, savename = 'Pattern_scaling'):

        fig, ax = plt.subplots(len(self.y_axis_variables), 1, figsize = (5, 5*len(self.y_axis_variables)))
        # Iterate over the y-axis variables
        for varind, y_var in enumerate(self.y_axis_variables):
            separate_pattern_flag = 0
            y_axis = []
    
        # Iterate over the benchmarking results of each target pattern
            for pattern_results in benchmarking_results:
                
                # If the y-axis variable is a list (e.g. filling fraction), take its average
                if type(pattern_results[y_var]) is list:
                    pattern_results[y_var] = np.mean(pattern_results[y_var])

                y_axis.append(pattern_results[y_var])
                separate_pattern_flag += 1

                # If all the results of the algorithm are collected, plot the results
                if separate_pattern_flag % len(x_axis) == 0:
                    try:
                        ax[varind].scatter(x_axis, y_axis, marker = 'o', label = CONFIGURATION_PLOT_LABELS[pattern_results["target"]])
                    except TypeError:
                        ax.scatter(x_axis, y_axis, marker = 'o', label = CONFIGURATION_PLOT_LABELS[pattern_results["target"]])
                    y_axis = []
            try:
                ax[varind].set_xlabel(x_label)
                ax[varind].set_ylabel(y_var.capitalize())
                ax[varind].set_title(f"{title} - {y_var.capitalize()}")
                ax[varind].legend(loc='best')
            except TypeError:
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_var.capitalize())
                ax.set_title(f"{title} - {y_var.capitalize()}")
                ax.legend(loc='best')

        if save:
            plt.savefig(f'./figs/{savename}')


# Set up the algorithms, target configurations, and system sizes
class Benchmarking():
    """
    An environment for studying the performance of rearrangement algorithms.
    
    Can be used to compare the scaling behavior of different algorithms, compare the time it takes for a single algorithm to prepare different target configurations, etc.
    
    ## Parameters
    - `algos` (list of `Algorithm` objects):
        the algorithms to compare.
    - `figure_output` (`BenchmarkingFigure`):
        an object for plotting.
    - `target_configs` (list of `Configurations` objects): 
        the target configurations to prepare.
    - `sys_sizes` (range):
        lengths of the square arrays that you want to look at (sqrt(N), where N is the number of tweezer sites).
    - `exp_params` (`PhysicalParams`):
        error and experimental parameters.
    - `n_shots` (int):
        number of repetitions per (algorithm or target config) per system size.

    ## Example Usage
    
    Runs the benchmarking round and generates a figure.
        `instance = Benchmarking()`
        `instance.run()`
        `instance.show()`
    """
    
    def __init__(self, 
                 algos: list = [Algorithm()], 
                 target_configs: list = [Configurations.MIDDLE_FILL], 
                 error_models_list: list = [ZeroNoise()],
                 phys_params_list: list = [PhysicalParams()],
                 sys_sizes: list = list(range(10,16)), 
                 rounds_list: list = [1],
                 figure_output: BenchmarkingFigure = BenchmarkingFigure(),
                 n_shots: int = 100,
                 n_species: int = 1):
        # initializing the sweep modules
        self.algos, self.n_algos = algos, len(algos)
        self.target_configs, self.n_targets = target_configs, len(target_configs)
        self.system_size_range, self.n_sizes = sys_sizes, len(sys_sizes)
        self.error_models_list, self.n_models = error_models_list, len(error_models_list)
        self.phys_params_list, self.n_parsets = phys_params_list, len(phys_params_list)
        self.rounds_list, self.n_rounds = rounds_list, len(rounds_list)

        # initializing other variables
        self.n_shots = n_shots
        self.figure_output = figure_output
        self.tweezer_array = AtomArray(n_species = n_species)

    def _input_error_check(self):
        if self.figure_output.figure_type not in ["scale", "pattern", "hist"]:
                raise Exception("Invalid figure type. Choose from 'scale', 'pattern', 'hist'.")
        if len(self.algos) > 1:
            if (len(self.target_config) > 1 or len(self.exp_params_list) > 1):
                raise Exception(f"Too many benchmarking variables. There are {len(self.algos)} \
                                algorithms, {len(self.target_config)} desired configurations, \
                                and {len(self.exp_params_list)} sets of parameters.")
            elif self.figure_output.figure_type == "pattern":
                raise Exception("Too many algorithms for scaling benchmarking of different target configurations.")
        else:
            if self.figure_output.figure_type == "scale" and \
               (len(self.target_configs) > 1 or len(self.exp_params_list) > 1):
                raise Exception("Cannot have more than 1 target config AND more than  for scaling benchmarking of different algorithms.")
            
            if self.figure_output.figure_type == "pattern" and len(self.algos) > 1:
                raise Exception("Too many algorithms for scaling benchmarking of different target configurations.")
            
            if self.figure_output.figure_type == "hist" and len(self.target_config) > 1:
                raise Exception("Too many target configurations for histogram benchmarking of different algorithms.")

    def plot_results(self, save = False, savename = None):
        if self.figure_output.figure_type == "scale":
            if savename == None:
                savename = 'scaling'
            self.figure_output.generate_scaling_figure(list(self.system_size_range), self.benchmarking_results, "Benchmarking results", "Array length (# atoms)", savename=savename, save=save)

        elif self.figure_output.figure_type == "hist":
            if savename == None:
                savename = 'histogram'
            self.figure_output.generate_histogram_figure(self.benchmarking_results, "Benchmarking results", "Array length (# atoms)")

        elif self.figure_output.figure_type == "pattern":
            if savename == None:
                savename = 'pattern'
            self.figure_output.generate_pattern_figure(list(self.system_size_range), self.benchmarking_results, "Benchmarking results", "Array length (# atoms)")


    def save(self, savename):
        if savename[-3:] == '.nc':
            savename = savename[0:-3]
        self.benchmarking_results.to_netcdf(f'data/{savename}.nc')
        print(f'Benchmarking object saved to `data/{savename}.nc`')
    
    def load(self, loadname):
        if loadname[-3:] == '.nc':
            loadname = loadname[0:-3]
        self.benchmarking_results = xr.open_dataset(f'data/{loadname}.nc', engine="netcdf4")
        print(f'Data from `data/{loadname}.nc` loaded to `self.benchmarking_results`.')

    def load_params_from_dataset(self, dataset: xr.Dataset):
        """ 
        Overwrites current parameters for benchmarking sweeps with those 
        from another xarray.Dataset object (e.g. `self.benchmarking_results`)
        
        Useful when wanting to retake data or play around with slightly different parameters.
        Also useful in recreating the figures from the atommovr paper.
        """
        self.algos = dataset['algorithm'].values
        self.target_configs = dataset['target'].values
        self.system_size_range = dataset['sys size'].values
        self.error_models_list = dataset['error model'].values
        self.phys_params_list = dataset['physical params'].values
        rounds_list = dataset['num rounds'].values
        self.rounds_list = []
        for round in rounds_list:
            self.rounds_list.append(int(round))
        self.n_shots = len(dataset['filling fraction'].values[0][0][0][0][0][0])

    def set_observables(self, observables: list):
        self.figure_output.y_axis_variables = observables


    def get_result_array_dims(self):
        """
        Updates the size and shape of the storage array
        based on the current set of parameters.
        """
        self.n_algos = len(self.algos)
        self.n_targets = len(self.target_configs)
        self.n_sizes = len(self.system_size_range)
        self.n_models = len(self.error_models_list)
        self.n_parsets = len(self.phys_params_list)
        self.n_rounds = len(self.rounds_list)

    def run(self, do_ejection: bool = False):
        """
        Run a round of benchmarking according to the parameters passed to the `Benchmarking()` object. 
        
        Saves the results in the variable `self.benchmarking_results`.
        """
        # Check the input conditions and raise an error if they are not met
        # self._input_error_check()

        # initializing result arrays
        self.get_result_array_dims()
        result_array_dims = [self.n_algos, 
                             self.n_targets, 
                             self.n_sizes, 
                             self.n_models, 
                             self.n_parsets, 
                             self.n_rounds]
        success_rate_array = np.zeros(result_array_dims, dtype = 'float')
        time_array = np.zeros(result_array_dims, dtype = 'float')
        fill_fracs_array = np.zeros(result_array_dims, dtype = 'object')
        wrong_places_array = np.zeros(result_array_dims, dtype = 'object')
        n_atoms_array = np.zeros(result_array_dims, dtype = 'object')
        n_targets_array = np.zeros(result_array_dims, dtype = 'object')

        # for xarray object
        dims = ("algorithm", "target", "sys size", "error model", "physical params", "num rounds")
        coords = {"algorithm": self.algos, 
                  "target": self.target_configs,
                  "sys size": self.system_size_range,
                  "error model": self.error_models_list,
                  "physical params": self.phys_params_list,
                  "num rounds": self.rounds_list}
        
        # iterating through sweep parameters and running benchmarking rounds
        for param_ind, parset in enumerate(self.phys_params_list):
            self.tweezer_array.params = parset
        
            self.init_config_storage = generate_random_init_configs(self.n_shots, 
                                                                    load_prob = self.tweezer_array.params.loading_prob, 
                                                                    max_sys_size = np.max(self.system_size_range),
                                                                    n_species = self.tweezer_array.n_species)
            for targ_ind, target in enumerate(self.target_configs):
                if target == Configurations.RANDOM:
                    self.target_config_storage = generate_random_target_configs(self.n_shots,
                                                                                load_prob = self.tweezer_array.params.target_occup_prob)
                for model_ind, error_model in enumerate(self.error_models_list):
                    self.tweezer_array.error_model = error_model

                    for size_ind, size in enumerate(self.system_size_range):
                        self.tweezer_array.shape = [size,size]

                        for alg_ind, algo in enumerate(self.algos):

                            for round_ind, num_rounds in enumerate(self.rounds_list):
                                success_rate, mean_success_time, fill_fracs, wrong_places, atoms_in_arrays, atoms_in_target = self._run_benchmark_round(algo, do_ejection=do_ejection, pattern = target, num_rounds=num_rounds)
                                # populating result arrays
                                success_rate_array[alg_ind, targ_ind, size_ind, model_ind, param_ind, round_ind] = success_rate
                                time_array[alg_ind, targ_ind, size_ind, model_ind, param_ind, round_ind] = mean_success_time
                                fill_fracs_array[alg_ind, targ_ind, size_ind, model_ind, param_ind, round_ind] = fill_fracs
                                wrong_places_array[alg_ind, targ_ind, size_ind, model_ind, param_ind, round_ind] = wrong_places
                                n_atoms_array[alg_ind, targ_ind, size_ind, model_ind, param_ind, round_ind] = atoms_in_arrays
                                n_targets_array[alg_ind, targ_ind, size_ind, model_ind, param_ind, round_ind] = atoms_in_target
                        # print(f'size {size} completed')
        
        success_rates_da = xr.DataArray(success_rate_array, dims=dims, coords = coords)
        success_times_da = xr.DataArray(time_array, dims=dims, coords = coords)
        fill_fracs_da = xr.DataArray(fill_fracs_array, dims=dims, coords = coords)
        wrong_places_da = xr.DataArray(wrong_places_array, dims=dims, coords = coords)
        n_atoms_da = xr.DataArray(n_atoms_array, dims=dims, coords = coords)
        n_targets_da = xr.DataArray(n_targets_array, dims=dims, coords = coords)
        
        self.benchmarking_results = xr.Dataset({'success rate': success_rates_da, 
                                                'time': success_times_da, 
                                                'filling fraction': fill_fracs_da,
                                                'wrong places': wrong_places_da,
                                                'n atoms': n_atoms_da,
                                                'n targets': n_targets_da})

            
    def _run_benchmark_round(self, algorithm, do_ejection: bool = False, pattern = None, num_rounds = 1) -> tuple[float, float, list, list, list]:
        success_times = []
        success_flags = []
        filling_fractions = []
        wrong_places = []
        atoms_in_arrays = []
        atoms_in_targets = []

        if pattern != Configurations.RANDOM:
            self.tweezer_array.generate_target(pattern, occupation_prob = self.tweezer_array.params.loading_prob)
                
        for shot in range(self.n_shots):
            # getting initial and final target configs
            initial_config = self.init_config_storage[shot][:self.tweezer_array.shape[0], :self.tweezer_array.shape[1]]
            self.tweezer_array.matrix = initial_config.reshape([self.tweezer_array.shape[0], self.tweezer_array.shape[1], self.tweezer_array.n_species])
            if pattern == Configurations.RANDOM:
                self.tweezer_array.target = self.target_config_storage[shot][:self.tweezer_array.shape[0], :self.tweezer_array.shape[1]].reshape([self.tweezer_array.shape[0], self.tweezer_array.shape[1], 1])
            init_count = 0
            while np.sum(initial_config) < np.sum(self.tweezer_array.target) and init_count < 100:
                self.tweezer_array.load_tweezers()
                initial_config = self.tweezer_array.matrix
                init_count += 1
            if init_count == 100:
                print(f'WARNING: could not find initial configuration with enough atoms ({np.sum(self.tweezer_array.target)}) in target). Consider aborting run and choosing more suitable parameters.')
            
            round_count = 0
            if num_rounds == 0 or not isinstance(num_rounds, int):
                raise ValueError(f'Number of rearrangement rounds (entered as {num_rounds}) cannot be 0 or a non-integer value.')
            while round_count < num_rounds:
                # generating and evaluating moves
                if self.tweezer_array.n_species == 1:
                    _, move_list, algo_success_flag = algorithm.get_moves(self.tweezer_array, do_ejection = do_ejection)
                else:
                    _, move_list, algo_success_flag = algorithm.get_moves(self.tweezer_array)
                t_total, _ = self.tweezer_array.evaluate_moves(move_list)
                success_flag = Algorithm.get_success_flag(self.tweezer_array.matrix,self.tweezer_array.target, do_ejection=do_ejection, n_species=self.tweezer_array.n_species)
                if success_flag == 1:
                    break
                round_count += 1

            success_flags.append(success_flag)
            if success_flag:
                success_times.append(t_total)

            # calculate filling fraction
            filling_fraction_config = np.multiply(self.tweezer_array.matrix, self.tweezer_array.target)
            filling_fractions.append(float(np.sum(filling_fraction_config)/np.sum(self.tweezer_array.target)))

            # Identify wrong places (atoms that are not in the target configuration)
            if do_ejection:
                wrong_places.append(int(np.sum(np.abs(self.tweezer_array.matrix - self.tweezer_array.target))))
            else:
                start_row, end_row, start_col, end_col = get_effective_target_grid(self.tweezer_array.target)
                wrong_places.append(int(np.sum(np.abs(self.tweezer_array.matrix[start_row:end_row, start_col:end_col] - self.tweezer_array.target[start_row:end_row, start_col:end_col]))))
            # Count atoms in array
            atoms_in_arrays.append(int(np.sum(self.tweezer_array.matrix)))
            atoms_in_targets.append(int(np.sum(self.tweezer_array.target)))

        return float(np.mean(success_flags)), float(np.mean(success_times)), filling_fractions, wrong_places, atoms_in_arrays, atoms_in_targets
