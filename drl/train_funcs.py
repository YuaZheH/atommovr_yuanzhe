import optuna
import gymnasium as gym
from stable_baselines3.common.callbacks import EvalCallback
from typing import Any, Dict, Optional

from atommover.utils.move_utils import get_move_list_from_AOD_actions
from atommover.utils.animation import single_species_image
from drl.policy_nns import CustomCNN

def visualize_testing_episode(env, model, use_action_masks: bool = True):  
    obs = env.reset()[0]
    count = 0
    dones = False
    print('Target:')
    single_species_image(env.target.reshape(env.array_len_y, env.array_len_x), padding=[[env.padding_y0, env.padding_y1], [env.padding_x0, env.padding_x1]])
    print('t=0')
    single_species_image(obs[0,:,:].reshape(env.array_len_y, env.array_len_x), padding=[[env.padding_y0, env.padding_y1], [env.padding_x0, env.padding_x1]])
    while not dones:
        obs = obs.reshape(2, env.array_len_y, env.array_len_x)
        if use_action_masks:
            action, _states = model.predict(obs, action_masks = env.get_action_mask())
        else:
            action, _states = model.predict(obs)
        obs, rewards, dones, truncs, info = env.step(action)
        horizontal_AOD_cmds = action[:env.array_len_x]
        vertical_AOD_cmds = action[env.array_len_x:]
        move_list = get_move_list_from_AOD_actions(horizontal_AOD_cmds, vertical_AOD_cmds)
        count += 1
        print(f't={count}, reward={rewards}')
        single_species_image(obs[0,:,:].reshape(env.array_len_y, env.array_len_x), move_list, padding=[[env.padding_y0, env.padding_y1], [env.padding_x0, env.padding_x1]])


## Hyperparameter Tuning with Optuna

def sample_ppo_params(trial: optuna.Trial) -> Dict[str, Any]:
    """
    Sampler for MaskablePPO hyperparams.

    :param trial:
    :return:
    """
    batch_size = trial.suggest_categorical("batch_size", [8, 16, 32, 64, 128, 256, 512])
    n_steps = trial.suggest_categorical("n_steps", [8, 16, 32, 64, 128, 256, 512, 1024, 2048])
    gamma = 1
    learning_rate = trial.suggest_float("learning_rate", 1e-5, 1, log=True)
    ent_coef = trial.suggest_float("ent_coef", 0.00000001, 0.1, log=True)
    clip_range = trial.suggest_categorical("clip_range", [0.1, 0.2, 0.3, 0.4])
    n_epochs = trial.suggest_categorical("n_epochs", [1, 5, 10, 20])
    gae_lambda = trial.suggest_categorical("gae_lambda", [0.8, 0.9, 0.92, 0.95, 0.98, 0.99, 1.0])
    max_grad_norm = trial.suggest_categorical("max_grad_norm", [0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 5])
    vf_coef = trial.suggest_float("vf_coef", 0, 1)
    # net_arch_type = trial.suggest_categorical("net_arch", ["tiny", "small", "medium"])
    # Uncomment for gSDE (continuous actions)
    # log_std_init = trial.suggest_float("log_std_init", -4, 1)
    # Uncomment for gSDE (continuous action)
    # sde_sample_freq = trial.suggest_categorical("sde_sample_freq", [-1, 8, 16, 32, 64, 128, 256])
    # Orthogonal initialization
    ortho_init = False
    normalize_images = False
    # ortho_init = trial.suggest_categorical('ortho_init', [False, True])
    # activation_fn_name = trial.suggest_categorical('activation_fn', ['tanh', 'relu', 'elu', 'leaky_relu'])
    # activation_fn_name = trial.suggest_categorical("activation_fn", ["tanh", "relu"])
    # lr_schedule = "constant"
    # Uncomment to enable learning rate schedule
    # lr_schedule = trial.suggest_categorical('lr_schedule', ['linear', 'constant'])
    # if lr_schedule == "linear":
    #     learning_rate = linear_schedule(learning_rate)

    # TODO: account when using multiple envs
    if batch_size > n_steps:
        batch_size = n_steps

    # Independent networks usually work best
    # when not working with images
    # net_arch = {
    #     "tiny": dict(pi=[64], vf=[64]),
    #     "small": dict(pi=[64, 64], vf=[64, 64]),
    #     "medium": dict(pi=[256, 256], vf=[256, 256]),
    #     "long": dict(pi=[64, 64, 64], vf=[64, 64, 64]),
    # }[net_arch_type]

    # activation_fn = {"tanh": nn.Tanh, "relu": nn.ReLU, "elu": nn.ELU, "leaky_relu": nn.LeakyReLU}[activation_fn_name]

    return {
        "n_steps": n_steps,
        "batch_size": batch_size,
        "gamma": gamma,
        "learning_rate": learning_rate,
        "ent_coef": ent_coef,
        "clip_range": clip_range,
        "n_epochs": n_epochs,
        "gae_lambda": gae_lambda,
        "max_grad_norm": max_grad_norm,
        "vf_coef": vf_coef,
        # "sde_sample_freq": sde_sample_freq,
        "policy_kwargs": dict(
            # log_std_init=log_std_init,
            # net_arch=net_arch,
            ortho_init=ortho_init,
            normalize_images = normalize_images,
            features_extractor_class=CustomCNN,
            features_extractor_kwargs=dict(features_dim=128),
        ),
    }



class TrialEvalCallback(EvalCallback):
    """
    Callback used for evaluating and reporting a trial.
    """

    def __init__(
        self,
        eval_env: gym.Env,
        trial: optuna.Trial,
        n_eval_episodes: int = 5,
        eval_freq: int = 10000,
        deterministic: bool = True,
        verbose: int = 0,
        best_model_save_path: Optional[str] = None,
        log_path: Optional[str] = None,
    ) -> None:
        super().__init__(
            eval_env=eval_env,
            n_eval_episodes=n_eval_episodes,
            eval_freq=eval_freq,
            deterministic=deterministic,
            verbose=verbose,
            best_model_save_path=best_model_save_path,
            log_path=log_path,
        )
        self.trial = trial
        self.eval_idx = 0
        self.is_pruned = False

    def _on_step(self) -> bool:
        if self.eval_freq > 0 and self.n_calls % self.eval_freq == 0:
            # Evaluate policy (done in parent class)
            super()._on_step()
            self.eval_idx += 1
            # report best or report current ?
            # report num_timesteps or elasped time ?
            self.trial.report(self.last_mean_reward, self.eval_idx)
            # Prune trial if need
            if self.trial.should_prune():
                self.is_pruned = True
                return False
        return True