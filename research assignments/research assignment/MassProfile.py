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
        self.filename = "%s_"%(galaxy)+ilbl+'.txt'
     
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
        print(COM_p)
        
        
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
        
        print(Vcirc)
        
        return Vcirc
    


# # Milky Way Mass Profile

# In[72]:


MW = MassProfile("MW", 0) # initializing MW classs


# In[73]:


r = np.arange(0.1, 30.5, 0.5) # radius range we are using up to 30kpc


# calculating enclosed mass for each galaxy component for each radius
enclosed_mass_MW_halo = MW.MassEnclosed(1,r)
enclosed_mass_MW_disk = MW.MassEnclosed(2,r)
enclosed_mass_MW_bulge = MW.MassEnclosed(3,r)

# calculating the total mass enclosed at each radius
total_mass_MW = MW.MassEnclosedTotal(r)


# In[157]:


a = 40 # the estimated scale factor
hernquist_profile = MW.HernquistMass(r, a, sum(MW.m_new*1e12)) # finding the hernquist model halo mass enclosed


# In[158]:


# plotting MW mass profile
plt.semilogy(r, total_mass_MW, label='total mass of MW', c='teal')
plt.semilogy(r, enclosed_mass_MW_halo, label='MW halo', c='midnightblue')
plt.semilogy(r, enclosed_mass_MW_disk, label='MW disk', c='m')
plt.semilogy(r, enclosed_mass_MW_bulge, label='MW bulge', c='sienna')
plt.semilogy(r,hernquist_profile, '--',label='Hernquist model', c='rosybrown', linewidth=3)

plt.ylabel('log(Mass in $M_{sol}$)')
plt.xlabel('radius in kpc')
plt.title("Milky Way Mass Profile at a = %s" % a)
plt.legend()


# In[62]:


G


# In[63]:


# calculating the circular velocity at each radius for each galaxy component
velocity_MW_halo = MW.CircularVelocity(1, r, G)
velocity_MW_disk = MW.CircularVelocity(2, r, G)
velocity_MW_bulge = MW.CircularVelocity(3, r, G)

# calculating the circular velocity at each radius for total mass
total_velocity_MW = MW.CircularVelocityTotal(r, G)


# # Milky Way Rotation Curve

# In[155]:


a = 40 # scale factor for MW

# finding the hernquist model of circular velocity
hernquist_profile = MW.HernquistVcirc(r,a, sum(MW.m_new*1e12), G) 


# In[156]:


# plotting MW rotation curves
plt.plot(r, total_velocity_MW, label='total mass of MW', c='teal')
plt.plot(r, velocity_MW_halo, label='MW halo', c='midnightblue')
plt.plot(r, velocity_MW_disk, label='MW disk', c='m')
plt.plot(r, velocity_MW_bulge, label='MW bulge', c='sienna')
plt.plot(r,hernquist_profile, '--',label='Hernquist Model', c='rosybrown', linewidth=3)
plt.yscale('log')
plt.ylabel('circular velocity in km/s')
plt.xlabel('radius in kpc')
plt.title("Milky Way Rotation Curve at a = %s" % a)
plt.legend()


# In[40]:


# initializing the M31 classs
M31 = MassProfile("M31", 0)
r = np.arange(0.25, 30.5, 0.5);print(r)

# calculating the enclosed mass for the galaxies component, at each radius
enclosed_mass_M31_halo = M31.MassEnclosed(1,r)
enclosed_mass_M31_disk = M31.MassEnclosed(2,r)
enclosed_mass_M31_bulge = M31.MassEnclosed(3,r)

#calculating total enclosed mass at each radius
total_mass_M31 = M31.MassEnclosedTotal(r)


# # M31 Mass Profile

# In[160]:


a = 60 # estimated scale factor for M31
hernquist_profile_M31 = M31.HernquistMass(r, a, sum(M31.m_new*1e12)) # finding hernquist model halo enclosed mass


# In[161]:


# plotting M31 mass profile
plt.semilogy(r, total_mass_M31, label='total mass of M31', c='teal')
plt.semilogy(r, enclosed_mass_M31_halo, label='halo', c='midnightblue')
plt.semilogy(r, enclosed_mass_M31_disk, label='disk', c='m')
plt.semilogy(r, enclosed_mass_M31_bulge, label='bulge', c='sienna')
plt.semilogy(r,hernquist_profile_M31, '--',label='Hernquist Model', c='rosybrown', linewidth=3)

plt.ylabel('Mass in $M_{sol}$')
plt.xlabel('radius in kpc')
plt.title("M31 Mass Profile at a = %s" % a)
plt.legend()


# # M31 Rotation Curve

# In[162]:


# calculating the circular velocity at each radius for each galaxy component
velocity_M31_halo = M31.CircularVelocity(1, r, G)
velocity_M31_disk = M31.CircularVelocity(2, r, G)
velocity_M31_bulge = M31.CircularVelocity(3, r, G)
# calculating the circular velocity at each radius for total mass
total_velocity_M31 = M31.CircularVelocityTotal(r, G)


# In[163]:


a = 60 #scale factor for M31
# finding the hernquist model of circular velocity
hernquist_profile_M31 = M31.HernquistVcirc(r,a, sum(M31.m_new*1e12), G)


# In[164]:


# plotting M31 rotation curves
plt.plot(r, total_velocity_M31, label='total mass of M31', c='teal')
plt.plot(r, velocity_M31_halo, label='M31 halo', c='midnightblue')
plt.plot(r, velocity_M31_disk, label='M31 disk', c='m')
plt.plot(r, velocity_M31_bulge, label='M31 bulge', c='sienna')
plt.plot(r,hernquist_profile_M31, '--',label='Hernquist Model', c='rosybrown', linewidth=3)
plt.yscale('log')
plt.ylabel('circular velocity in km/s')
plt.xlabel('radius in kpc')
plt.title("M31 rotation curve at a = %s" % a)
plt.legend()


# # M33 Mass Profile

# In[46]:


# initializing the M33 classs
M33 = MassProfile("M33", 0)
r = np.arange(0.25, 30.5, 0.5);print(r)

# calculating enclosed mass for each galaxy component for each particular radius
enclosed_mass_M33_halo = M33.MassEnclosed(1,r)
enclosed_mass_M33_disk = M33.MassEnclosed(2,r)
#calculating the total mass enclosed at each radius
total_mass_M33 = M33.MassEnclosedTotal(r)


# In[165]:


a = 80 # estimated scale factor
hernquist_profile_M33 = M33.HernquistMass(r, a, sum(M33.m_new*1e12)) # finding hernquist model halo mass encloseda = 80 # estimated scale factor
hernquist_profile_M33 = M33.HernquistMass(r, a, sum(M33.m_new*1e12)) # finding hernquist model halo mass enclosed


# In[166]:


# plotting the M33 mass profile
plt.semilogy(r, total_mass_M33, label='total mass of M33', c='teal')
plt.semilogy(r, enclosed_mass_M33_halo, label='halo', c='midnightblue')
plt.semilogy(r, enclosed_mass_M33_disk, label='disk', c='m')
plt.semilogy(r,hernquist_profile_M33, '--',label='Hernquist Model', c='rosybrown', linewidth=3)

plt.ylabel('Mass in $M_{sol}$')
plt.xlabel('radius in kpc')
plt.title("M33 Mass Profile at a = %s" % a)
plt.legend()


# # M33 Rotation Curve

# In[49]:


# calculating the circular velocity at each particular radius for each galaxy component
velocity_M33_halo = M33.CircularVelocity(1, r, G)
velocity_M33_disk = M33.CircularVelocity(2, r, G)
# calculating the circular velocity at each radius for total mass
total_velocity_M33 = M33.CircularVelocityTotal(r, G)


# In[168]:


a = 80 # the scaled factor
# finding the hernquist model of circular velocity
hernquist_profile_M33 = M33.HernquistVcirc(r,a, sum(M33.m_new*1e12), G)


# In[170]:


# plotting the M33 rotation curve
plt.plot(r, total_velocity_M33, label='total mass of M33', c='teal')
plt.plot(r, velocity_M33_halo, label='M33 halo', c='midnightblue')
plt.plot(r, velocity_M33_disk, label='M33 disk', c='m')
plt.plot(r,hernquist_profile_M33, '--',label='Hernquist Model', c='rosybrown', linewidth=3)
plt.yscale('log')
plt.ylabel('circular velocity in km/s')
plt.xlabel('radius in kpc')
plt.title("M33 Rotation Curve at a = %s" % a)
plt.legend()


# In[ ]:




