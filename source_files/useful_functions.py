import numpy as np
import simulationParameters as para

quadrat_lenght = para.lenght_of_a_quadrat*1000 #quadrat lenght in m

x_min = 55000-quadrat_lenght*para.nb_fic_x_left
y_min = 1708000-quadrat_lenght*para.nb_fic_y_bottom

def transform_data(data, index_x, index_y):
    '''Takes an array data such that (data[i,index_x], data[i,index_y])
    is the extended lambert II coordinates of the i-th data and returns the 
    same data where the extended lambert II coordinates have been transformed 
    in quadrat indexes. Moreover, it excludes the data of Corsica.'''
    z = data[data[:,index_x]<1031000+16000,:] #exclude Corsica
    z[:,index_x] = np.round((z[:,index_x]-x_min)/quadrat_lenght)
    z[:,index_y] = np.round((z[:,index_y]-y_min)/quadrat_lenght)
    return z

