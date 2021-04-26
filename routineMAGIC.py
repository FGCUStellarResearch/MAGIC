from functionMAGICfieldCR import MAGIC_Field_CR
from astropy.table import QTable, Table, Column
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import astropy.units as u
for ra_pos in np.linspace(0,0.1,2):
    for dec_pos in np.linspace(0,0.1,2):
        MAGIC_Field_CR(ra_pos,dec_pos)
        
        
    
        

     
    