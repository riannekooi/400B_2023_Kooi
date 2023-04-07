#!/usr/bin/env python
# coding: utf-8

# In[12]:


import numpy as np
import astropy.units as u
# import Latex module so we can display the results with symbols
from IPython.display import Latex
from IPython.display import display 


# In[13]:


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


# In[14]:


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


# In[15]:


# MW: Compute Mass for each component
MW_halo = ComponentMass("MW_000.txt",1)
MW_disk = ComponentMass("MW_000.txt",2)
MW_bulge = ComponentMass("MW_000.txt",3)


# In[16]:


display(Latex(
    r"MW halo Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(MW_halo)))
display(Latex(
    r"MW disk Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(MW_disk)))
display(Latex(
    r"MW bulge Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(MW_bulge)))


# In[17]:


# Total MW Mass 
MW_total = MW_halo + MW_disk + MW_bulge
# MW Baryon Fraction
MW_f_bar = (MW_disk + MW_bulge) / MW_total


# In[18]:


display(Latex(
    r"MW total Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(MW_total)))
display(Latex(
    r"MW baryon fraction $f_{{bar}} = {:.3f}$".format(MW_f_bar)))


# In[20]:


# M31: Compute Mass for each component
M31_halo = ComponentMass("M31_000.txt",1)
M31_disk = ComponentMass("M31_000.txt",2)
M31_bulge = ComponentMass("M31_000.txt",3)


# In[21]:


display(Latex(
    r"M31 Halo Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(M31_halo)))
display(Latex(
    r"M31 Disk Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(M31_disk)))
display(Latex(
    r"M31 Bulge Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(M31_bulge)))


# In[22]:


# Total M31 Mass
M31_total = M31_halo + M31_disk + M31_bulge
# M31 Baryon Fraction 
M31_f_bar = (M31_disk + M31_bulge) / M31_total


# In[23]:


display(Latex(
    r"M31 Total Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(M31_total)))
display(Latex(
    r"M31 baryon fraction $f_{{bar}} = {:.3f}$".format(M31_f_bar)))


# In[24]:


# M33: Compute Mass for each component
M33_halo = ComponentMass("M33_000.txt",1)
M33_disk = ComponentMass("M33_000.txt",2)


# In[25]:


display(Latex(
    r"M33 Halo Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(M33_halo)))
display(Latex(
    r"M33 Disk Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(M33_disk)))


# In[26]:


# Total M33 Mass
M33_total = M33_halo + M33_disk
# M33 Baryon Fraction 
M33_f_bar = M33_disk / M33_total


# In[27]:


display(Latex(
    r"M33 Total Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(M33_total)))
display(Latex(
    r"M33 baryon fraction $f_{{bar}} = {:.3f}$".format(M33_f_bar)))


# In[28]:


# Total mass for the Local Group
LG_total = MW_total + M31_total + M33_total
display(Latex(
    r"Local Group Mass: ${:.3f}\times 10^{{12}}M_{{\odot}}$".format(LG_total)))


# In[29]:


# Baryon fraction for the Local Group 
LG_f_bar = (MW_disk + MW_bulge + M31_disk + M31_bulge + M33_disk) / LG_total


# In[30]:


# print out table
print()
print("Galaxy Name  | Halo Mass   |  Disk Mass   | Bulge Mass  | Total Mass  | f_bar ")
print("             | [1e12 Msun] |  [1e12 Msun] | [1e12 Msun] | [1e12 Msun] |       ")
print("-------------|-------------|--------------|-------------|-------------|-------")
print(" Milky Way   | {:<8.3f}    | {:<8.3f}     | {:<8.3f}    | {:<8.3f}    | {:<8.3f}".format(MW_halo, MW_disk, MW_bulge, MW_total, MW_f_bar))
print(" M31         | {:<8.3f}    | {:<8.3f}     | {:<8.3f}    | {:<8.3f}    | {:<8.3f}".format(M31_halo, M31_disk, M31_bulge, M31_total, M31_f_bar))
print(" M33         | {:<8.3f}    | {:<8.3f}     | -           | {:<8.3f}    | {:<8.3f}".format(M33_halo, M33_disk, M33_total, M33_f_bar))
print(" Local Group | -           | -            | -           | {:<8.3f}    | {:<8.3f}".format(LG_total, LG_f_bar))
print()


# In[ ]:




