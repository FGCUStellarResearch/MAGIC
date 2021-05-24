from functionMAGICfieldCR import MAGIC_Field_CR
from astropy.table import QTable, Table, Column
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import astropy.units as u
import time

test = np.empty((8,))

for ra_index, ra_pos in enumerate(np.linspace(360,360,1)):
    for dec_index, dec_pos in enumerate(np.linspace(50,60,2)):
        time.sleep(5)
        #test[ra_index,dec_index] = MAGIC_Field_CR(ra_pos,dec_pos)
        print(ra_pos)
        print(dec_pos)
        new_add = MAGIC_Field_CR(ra_pos,dec_pos)
        #test = np.append(test,new_add,axis=0)
        test = np.vstack((test,new_add))

print(test)

        
    
        

     
    
