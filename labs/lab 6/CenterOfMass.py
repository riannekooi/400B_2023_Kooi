#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Homework 4
# Center of Mass Position and Velocity
# Rianne Kooi


# In[4]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[5]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        # write your own code to compute the generic COM 
        #using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        a_com = np.sum(a*m)/np.sum(m)
        # ycomponent Center of mass
        b_com = np.sum(b*m)/np.sum(m)
        # zcomponent Center of mass
        c_com = np.sum(c*m)/np.sum(m)
        
        # return the 3 components separately
        return a_com, b_com, c_com
    
    
    def COM_P(self, delta=0.1):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2)
        
        print('r_com', r_COM)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(r_new < r_max)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2, y2, z2, m2)
            # compute the new 3D COM position
            # write your own code below
            r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2)
            
            print('r com: ', r_COM2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)
            # uncomment the following line if you want to check this                                                                                               
            # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : r_max, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            r_max /= 2.0
            # check this.                                                                                              
            #print ("maxR", r_max)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            x_new = self.x - x_COM2
            y_new = self.y - y_COM2
            z_new = self.z - z_COM2
            r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)
            
            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # create an array (np.array) to store the COM position                                                                                                                                                       
            p_COM = np.array([x_COM, y_COM, z_COM])

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        self.r_com = np.sqrt(p_COM[0]**2 + p_COM[1]**2 + p_COM[2]**2)
        return np.round(p_COM*u.kpc, 2)
        
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.

        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # the max distance from the center that we will use to determine 
        #the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        # write your own code below
        xV = self.x - (x_COM/u.kpc)
        yV = self.y - (y_COM/u.kpc)
        zV = self.z - (z_COM/u.kpc)
        rV = np.sqrt(xV**2 + yV**2 + zV**2)*u.kpc
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(rV<rv_max)
        
        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vx_new = self.vx[indexV]
        vy_new = self.vy[indexV]
        vz_new = self.vz[indexV]
        m_new =  self.m[indexV]
        
        # compute the center of mass velocity using those particles
        # write your own code below
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)
        
        # create an array to store the COM velocity
        # write your own code below
        v_COM = np.array([vx_COM, vy_COM, vz_COM])

        # return the COM vector
        # set the correct units usint astropy
        # round all values
        return np.round(v_COM*u.km/u.s, 2)
     
    


# In[6]:


# Create a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MW_COM = CenterOfMass("MW_000.py", 2)


# In[7]:


# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM
MW_COM_p = MW_COM.COM_P(0.1)
print(MW_COM_p)
MW_COM_v = MW_COM.COM_V(MW_COM_p[0], MW_COM_p[1], MW_COM_p[2])
print(MW_COM_v)


# In[8]:


M31_COM = CenterOfMass("M31_000.py", 2)
M31_COM_p = M31_COM.COM_P(0.1)
print(M31_COM_p)
M31_COM_v = M31_COM.COM_V(M31_COM_p[0], M31_COM_p[1], M31_COM_p[2])
print(M31_COM_v)


# In[9]:


M33_COM = CenterOfMass("M33_000.py", 2)
M33_COM_p = M33_COM.COM_P(0.1)
print(M33_COM_p)
M33_COM_v = M33_COM.COM_V(M33_COM_p[0], M33_COM_p[1], M33_COM_p[2])
print(M33_COM_v)


# In[12]:


# Making a Table of Results

# we first gather all results and put them into a numpy array
# *tuple() will extract the tuple, kind of like C-language
tab_results = ['MW COM', *tuple(MW_COM_p.value), *tuple(MW_COM_v.value), 
               'M31 COM', *tuple(M31_COM_p.value), *tuple(M31_COM_v.value), 
               'M33 COM', *tuple(M33_COM_p.value), *tuple(M33_COM_v.value)]
tab_results = np.reshape(tab_results, (3, 7)) # convert 1D array to 2D

# make and display the table with astropy
t = tbl.Table(tab_results, 
              names = ['Galaxy', 'X-component (kpc)', 'Y-Component (kpc)', 
                       'Z-Component (kpc)','VX-component (km/s)', 
                       'VY-Component (km/s)', 'VZ-Component (km/s)'])
t.show_in_notebook()


# In[10]:


#What is the magnitude of the current separation (in kpc) and velocity (in km/s) between the MW and M31? 
#Round your answers to three decimal places. 

print('MW COM Position',MW_COM_p)
print('MW COM Position',M31_COM_p)

sep_p = np.sqrt((M31_COM_p[0]-MW_COM_p[0])**2 + (M31_COM_p[1]-MW_COM_p[1])**2 + (M31_COM_p[2]-MW_COM_p[2])**2)
print('separation:',np.round(sep_p, 3))

sep_v = np.sqrt((M31_COM_v[0]-MW_COM_v[0])**2 + (M31_COM_v[1]-MW_COM_v[1])**2 + (M31_COM_v[2]-MW_COM_v[2])**2)
print('velocity separation:',np.round(sep_v, 3))


# In[11]:


#What is the magnitude of the current separation (in kpc) and velocity (in km/s) between M33 and M31? 
#Round your answers to three decimal places.
print(M31_COM_p)
print(M33_COM_p)

sep_p = np.sqrt((M31_COM_p[0]-M33_COM_p[0])**2 + (M31_COM_p[1]-M33_COM_p[1])**2 + (M31_COM_p[2]-M33_COM_p[2])**2)
print('separation:',np.round(sep_p, 3))

sep_v = np.sqrt((M31_COM_v[0]-M33_COM_v[0])**2 + (M31_COM_v[1]-M33_COM_v[1])**2 + (M31_COM_v[2]-M33_COM_v[2])**2)
print('velocity separation:',np.round(sep_v, 3))


# In[ ]:


#Given that M31 and the MW are about to merge, why is the iterative process to determine the COM is important?

#The interative process to determine the COM is important because as time increases and M31 and the Milky Way 
# get closer and closer, then the Center of Mass will change extremely rapidly. Using the iterative process will 
# ensure that the true COM is calculated. 

