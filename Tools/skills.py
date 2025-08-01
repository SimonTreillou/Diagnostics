import numpy as np


def skill_Willmott(mod,obs):
    # For reference, see Rijnsdorp et al. 2017
    mean_obs = np.mean(obs)
    skill = 1 - np.sum((mod - obs)**2) / np.sum((np.abs(mod - mean_obs) + np.abs(obs - mean_obs))**2)
    return skill

def nrmse(mod,obs):
    # For reference, see Martins et al. 2022
    mean_obs = np.mean(obs)
    skill =100 * np.sqrt(np.nansum((obs-mod)**2)/np.nansum(obs**2))
    return skill