# -*- coding:Utf-8 -*-
"""
Enter in this script the name of the folders associated to computations of the
log-likelihoods on regulars grids that you want to plot in the variable 'simu':
    simu = ['name_folder1', 'name_folder2',...,'name_folderN'].
You can also change some graphical parameters (plot on log-scale, color, size of
labels,...). See parameters below.
Then run the script 'main_script.py' and select action 9.
The graphs will be plotted in the folder 'Simulations/Compare_MH_1' (or
in the folder 'Simulations/Compare_MH_(i+1)' if the folder 
'Simulations/Compare_MH_i' already exists).
"""

simus = []

mini = - 25150
non_test = mini
maxi = "None" #-24600

logscale_r_allee = True 
value_0 = 10**(-15)

plot_beta0 = ["all"]  
plot_D = ["all"]
plot_gamma = ["all"]
plot_kappa = ["all"]
plot_r_sym = ["all"]
plot_r_allee = ["all"]
plot_beta1 = ["all"]
plot_S_sat = ["all"]
plot_C_recover = ["all"]
plot_C_init = ["all"] 


fontsizex = '16' #fontsize of the xlabel
fontsizey = '16' #fontsize of the ylabel
fontsizecbar = '16' #fontsize of the label of the cbar
fontsizeticks = '10' #fontsize of the ticks
