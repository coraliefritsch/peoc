# -*- coding:Utf-8 -*-
"""
Enter in this script the parameters for plotting the evolution of the parameters
of different runs of the Metropolis-Hastings algorithm.
This action can be used for example to observe the 
convergence of parameters for the algorithm started from different initial 
conditions.
To plot the  the evolution of the parameters, add the names of the simulation 
folders in the variable 'simu':
    simu = ['name_folder1', 'name_folder2',...,'name_folderN'].
You can also change some graphical parameters (plot on log-scale, color, size of
labels,...). See parameters below.
Then run this script 'main_script.py' and select action 6.
The graphs will be plotted in the folder 'Simulations/Compare_MH_1' (or in the 
folder 'Simulations/Compare_MH_(i+1)' if the folder 
'Simulations/Compare_MH_i' already exists).
"""

simu = ['MH_initial_condition_1', 'MH_initial_condition_2',
        'MH_initial_condition_3', 'MH_initial_condition_4',
        'MH_initial_condition_5', 'MH_initial_condition_6',
        'MH_initial_condition_7', 'MH_initial_condition_8',
        'MH_initial_condition_9', 'MH_initial_condition_10']

format_ = ".png" #format of figures
evol2D = False #if True then plot the evolution of all pairs of variables

#if plot_log[i] = 1 then the i-th variable is plotted on a logarithmic scale
plot_log = [0,0,0,0,0,0,0,0,0,0]
#the i-th variable is plotted every k[i] iterations (to reduce size of files)
k = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

#size of the labels and ticks
fontsizelabelx = 14
fontsizelabely = 14
fontsizetickx = 10
fontsizeticky = 8

#color of the graph for each simulation
c = ['r', 'g', 'darkorange', 'k', 'c', 'm', 'y', 'b', 'indigo', 'gold',
         'aquamarine', 'lightblue', 'maroon', 'lightslategray',
         'darkslategray']
