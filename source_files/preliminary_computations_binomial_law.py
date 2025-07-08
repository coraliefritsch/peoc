# -*- coding:Utf-8 -*-
"""
Preliminary computations for the log likelihood calculus of the binomial laws:
This file compute clb_likeli defined by:
clb_likeli[a][:,0] = log[(observed trees number) choose (infected trees number)]
clb_likeli[a][:,1] = infected observed trees number of year a
clb_likeli[a][:,2] = non infected observed trees number of year a
"""
import numpy as np
import simulationParameters as para
import observations
from scipy.special import gammaln as gamln

clb_likeli = []
for year in observations.observation:
    obs_trees = para.number_of_observed_trees*year[3,:]#number of observed trees
    #number of infected and non infected observed trees:
    obs_infect_trees = np.int_(np.round(obs_trees*year[2,:]))
    obs_non_infect_trees = obs_trees - obs_infect_trees
    #binomial coefficients:
    combiln_obs = (gamln(obs_trees+1) - (gamln(obs_infect_trees+1)
                                            + gamln(obs_non_infect_trees+1)))
    z = np.zeros((np.size(year[3,:]),3))
    z[:,0] = combiln_obs
    z[:,1] = obs_infect_trees
    z[:,2] = obs_non_infect_trees

    clb_likeli.append(z)
    
del z, obs_trees, obs_infect_trees, obs_non_infect_trees, combiln_obs, year
