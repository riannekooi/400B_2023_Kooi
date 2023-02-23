#!/usr/bin/env python
# coding: utf-8

# In[56]:


#imports needed
import numpy as np
import astropy.units as u
import astropy.table as tbl
import matplotlib.pyplot as plt

from ReadFile import Read
from CenterOfMass import CenterOfMass
from astropy.constants import G


# In[71]:


class MassProfile:
# Class to define 

    def __init__(self, galaxy, snap):
        
        ''' 
        class to find the mass profile of a galaxy 
            
            parameters
            
                galaxy: 
                a string with galaxy name (like MW, M31, etc)
                snap: Snapshot number
            '''
        
        # adding a string of the filename to the value 000
        ilbl = '000'+str(snap)
        # removing all digits except the last 3 
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy)+ilbl+'.py'
     
        # reading the data in the given file using previously defined Read function
        self.time, self.total, self.data = Read(self.filename)                                                                                             

        # storing the positions, velocities, & mass
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.vx = self.data['vx']*u.km/u.s
        self.vy = self.data['vy']*u.km/u.s
        self.vz = self.data['vz']*u.km/u.s
        self.m = self.data['m']
        
        #set the galaxy name as a global property
        self.gname = galaxy
        
    def MassEnclosed(self, ptype, radii):
        '''
        this function will compute the enclosed mass 
        for a given radius of the COM position for a particular galaxy's component.

        inputs:
            ptype:
                the particle type to use for COM calculations (1, 2, 3)
            radii: 
                an array of radii in magnitude

        outputs: 
            masses_enclosed: 
                an array of enclosed mass at each particular radius specified
                in units of M_sun 

        '''
        #creating the CenterOfMass objects and calling previously defined COM_P
        COM = CenterOfMass(self.filename, ptype)
        COM_p = COM.COM_P(0.1)
        
        
        
        index = np.where(self.data['type'] == ptype)
        
        #loop over the radius array to define particles that are enclosed within the radius given
        
        self.m_new = self.m[index]
        
        masses_enclosed = np.zeros(len(radii))
        for i in range(len(radii)):
            
            distance_magnitude = np.sqrt((self.x[index]-COM_p[0])**2 
                                         +(self.y[index]-COM_p[1])**2 
                                         +(self.z[index]-COM_p[2])**2)
            
            enclosed_index = np.where(distance_magnitude/u.kpc < radii[i])
            
            #storing the sum of masses of the particles within the particular radius
            masses_enclosed[i]= sum(self.m_new[enclosed_index])
            
        return masses_enclosed * 1e10 *u.Msun
    
    def MassEnclosedTotal(self, radii):
        ''' this calculates the total enclosed mass from each particle type
        
        inputs:
            radii: 
                array of radii to calculate the enclosed mass at each point
        
        outputs:
            total_mass: 
                array of the toal mass enclosed at each radii in M_sun
        
        '''
        
        # calulating the mass enclosed for each particle type
        halo_mass = self.MassEnclosed(1, radii)
        disk_mass = self.MassEnclosed(2, radii)
        if self.gname != 'M33': # M33 doesn't have a buldge mass, so we skip it
            bulge_mass = self.MassEnclosed(3, radii)  
            total_mass = halo_mass + disk_mass + bulge_mass
        else:
            total_mass = halo_mass + disk_mass
            
        return total_mass
    
    def HernquistMass(self, radii, a, M_halo):
        '''
        computes the mass enclosed within a particular radius
        
        inputs:
            radius: 'array'
                an array of radii values
            a: 'float'
                hernquist scale factor
            M_halo: 'array'
                the halo mass of the galaxy
            
        outputs:
            M: the mass enclosed at a given radius in M_sun
            
        '''
        
        # hernquist mass profile
        
        M = (M_halo*radii**2/(a+radii)**2 )*u.Msun
            
        return M
    
    def CircularVelocity(self, ptype, radii, G):
        '''
        this function calculates the circular speed at each particular radii.
        
        inputs:
            ptype: 
                particle type (as a float)
            radii: 
                an array of radii for the circular speed to be calculated at 
            G: 
                gravitational constant in kms units
        
        output:
            Vcirc:
                an array of circular speeds with units km/s, rounded to two decimal places
        
        '''
        
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) # converting the units
        
        # finding mass enclosed
        M_enclosed = self.MassEnclosed(ptype, radii)
        
        # applying circular speed formula
        Vcirc = np.round(np.sqrt(G*M_enclosed/(radii*u.kpc)), 2)
        
        #print(Vcirc)
        
        return Vcirc
    
    def CircularVelocityTotal(self, radii, G):
        '''
        this function calculates the total circular velocity 
        representing the total Vcirc created by all the galaxy components
        at each particular radius from the input array
        
        inputs:
            radii:
                an array of radii
            
            G:
                gravitational constant in kms units
        
        outputs:
            total_Vcirc:
                an array of circular velocity (in units of km/s)
                
        '''
        
        # calculating circular velocity for each of the particle types
        halo_Vcirc = self.CircularVelocity(1, radii, G)
        disk_Vcirc = self.CircularVelocity(2, radii, G)
        if self.gname != 'M33': # M33 doesn't have a buldge mass, so we skip it
            bulge_Vcirc = self.CircularVelocity(3, radii, G)  
            total_Vcirc = halo_Vcirc + disk_Vcirc + bulge_Vcirc
        else:
            total_Vcirc = halo_Vcirc + disk_Vcirc
            
        return total_Vcirc
    
    def HernquistVcirc(self,radius,a, M_halo, G):
        '''
        this function computes the circular speed using the hernquist mass profile
        
        inputs:
            radius:
                an array of radii for the circular speed to be computed at
            a: 
                Hernquist scale factor (float)
            M_halo:
                an array of the halo mass of the galaxy
            
        
        '''
        
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) # converting units
        
        Mass_Hern = self.HernquistMass(radius,a,M_halo) # calculating the hernquist mass
        
        Vcirc = np.round(np.sqrt(G*Mass_Hern/(radius*u.kpc)), 2) #round and add units
        
        
        
        return Vcirc
    


