"""
Enter in this script the discretization of each variable for coputations of the
log-likelihood on a regular grid. Each set of values has to be an array. Then 
run this script 'main_script.py' and select action 8.
"""
import numpy as np

v_r_symp = np.linspace(1000, 1000, 1)
v_beta_1 = np.linspace(0.00 ,0.00 ,1)
v_S_sat = np.linspace(80.,100.,3)
v_C_recover = np.linspace(1.,1.,1)
v_C_init = np.linspace(0.008,0.009,3)
v_gamma = np.linspace(21.001,23.001,5)
v_beta_0 = np.linspace(2, 100, 50)
v_D = np.linspace(8.,35.,28)
v_kappa = np.linspace(0.025,0.25,10)
v_r_allee = np.linspace(0,0.0, 1)
#v_r_allee = np.append(np.linspace(0.,0.,1),10**(np.linspace(-10.,5.,16)))
    

num_workers = 2 #number of computation in parallel (number of cores)
