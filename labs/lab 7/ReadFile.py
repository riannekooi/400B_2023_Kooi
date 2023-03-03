#!/usr/bin/env python
# coding: utf-8

# In[20]:


#input: MW_000 file and read it
#output: the time and total number of particles as variables
#output: particle type, mass, x,y,z, vx, vy, vz as a data array


# In[21]:


import numpy as np
import astropy.units as u


# In[22]:


def Read(filename):
    file = open(filename, 'r')
    
#Row 1 is the time in units of Myr (SnapNumber*10/0.7)
    line1 = file.readline() #reading the first line of the file
    label, value = line1.split() #splitting the label of the line and actual values
    time = float(value)*u.Myr

    #Row 2 is the total number of particles
    line2 = file.readline()#reading the second line of the file
    label, value = line2.split() #splitting the label of the line and actual values
    numParticles = float(value)
    file.close()#closing the file that was being read
    
    #storing the rest of the data in the file
    #data stores the particle type, the mass, x,y,z coords, and the velocity in x, y, and z as a data array
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    return time, numParticles, data


# In[23]:


