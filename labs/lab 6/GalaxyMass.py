#!/usr/bin/env python
# coding: utf-8

# In[12]:


import numpy as np
import astropy.units as u
# import Latex module so we can display the results with symbols
from IPython.display import Latex
from IPython.display import display 
from ReadFile import Read


def ComponentMass(filename, part_type):
    """ Function to read the data from a given snapshot and return the total mass
    of the specified particle type.
    
    INPUTS
    ------
    filename: 'str'
        Name of the snapshot file to read
    part_type: 'int: 1,2,3'
        Particle type that will be summed to return mass
        
        
    OUTPUTS
    ------
    mass: 'float'
        Total mass of teh specified particle type in 1e12 solar masses
    """
  
    # read teh particle data from the specified file
    time, total, data = Read(filename)
    
    # select particles with the same type and sum up the mass
    mass = np.sum(data[data['type'] == part_type]['m'])
    
    # round and return the result in the correct units (1e12)
    return np.round(mass*1e10/1e12, 3)


