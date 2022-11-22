#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


def ET_PM(T,A_P,VPD,NETRAD,G,W):
    
    '''Input: T(Temperature[°C]), A_P (Atmospheric Pressure[kPa]), VPD(Vapour Pressure Deficit[kPa]), NETRAD (Net Radiation [MJm-2/day]),
    G( Ground Heat Flux [MJm-2/day]),W(Wind Speed [m/s]);
    Output: Evapotranspiration according to Penman-Monteith [mm/d'''
    
    ##Unit Conversions:
    ##W/m^2 = 0.0864 MJm-2/day
    ##J = 1Nm = 1kg m^2/s^2 = 1Pa * m^3
    ##Pa = 1N/m^2 = 1kg/ms^2= 1J / m^3
    ##1day = 86400s
    ##1s = 1 days/86400
    
    ET = np.zeros_like(T)
    for i in range(len(T)):
        slope = 4098*(0.6108*np.exp(17.27*T[i]/(T[i]+237.3))/(T[i]+237.3)**2)#[kPa/°C]
        c_p = 1.013 * 10**(-3)  #[MJ/kg*°C]
        gamma = 0.665*10**(-3)*A_P[i] #[kPa/°C]

        a = 0.408*slope*(NETRAD[i]-G[i]) #kPa/°C * MJ/m^2*day = kPa^2/°C^2 * m/day*1000
        b = gamma*900/(T[i]+273)*W[i]*VPD[i] #kPa^2/°C *m/s
        c = (slope+gamma*(1+0.34*W[i])) #kPa/°C
        
        # filters for index, that contains data
        if (T[i]!= -9999 and A_P[i]!= -9999 and P[i]!=-9999 and VPD[i]!=-999.9 and NETRAD[i]!=-863.9136000000001 and G[i] != -863.9136000000001 and W[i] != -863.9136000000001 ):
            ET[i] = (a+b)/c     
        else: 
            ET[i] = np.nan
            

    return ET


def ET_PT(T,A_P,NETRAD,G,alpha):
    
        
    '''Input: T(Temperature[°C]), A_P (Atmospheric Pressure[kPa]), NETRAD (Net Radiation [MJm-2/day]),
    G( Ground Heat Flux [MJm-2/day]), alpha(An empirical constant accounting for the vapor pressure deficit and resistance values);
    Output: Evapotranspiration according to Penman-Monteith in [mm/d] '''
    
    ET_PT = np.zeros_like(T)
    
    for i in range(len(T)):
        
        slope = 4098*(0.6108*np.exp(17.27*T[i]/(T[i]+237.3))/(T[i]+237.3)**2)#[kPa/°C]
        gamma = 0.665*10**(-3)*A_P[i] #[kPa/°C]
        
        # filters for index, that contains data
        if (T[i]!= -9999 and A_P[i]!= -9999  and NETRAD[i]!=-863.9136000000001 and G[i] != -863.9136000000001):
            ET_PT[i]= alpha* slope/(slope+gamma) *(NETRAD[i]-G[i])
        else: 
            ET_PT[i]= np.nan
        
    
    return ET_PT


# In[ ]:




