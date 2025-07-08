# -*- coding:Utf-8 -*-
"""
Enter in this script the parameters for plotting the posterior distributions of
different AMIS runs (action 7 of main_script). This action can be used for 
example to compare the posterior distributions for different climate variables.
To plot the posterior distributions, add 
- the names of the simulation folders in the variable 'simu': which has to be 
    simu = ['name_folder1', 'name_folder2',...,'name_folderN']
- the labels that have to be indicated on the graphs in the variable 'l':
       l = ['label_simu1','label_simu2',...,'label_simuN']
- the color of the plot for each simulation in the variable 'c' : 
       l = ['color_simu1','color_simu2',...,'color_simuN']
then run this script 'main_script.py' and select action 7.
The graphs will be plotted in the folder 'Simulations/Compare_posterior_1' (or
in the folder 'Simulations/Compare_posterior_(i+1)' if the folder 
'Simulations/Compare_posterior_i' already exists).
"""

#names of the simulation folders
simu = ['AMIS_with_temperature_index_T26', 'AMIS_with_temperature_index_T28']
#legend labels of the simulations
l = ['T26', 'T28']
#color of the simulations
c = ['b', 'r', 'g', 'c', 'y', 'm', 'k']
