#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas
import zipfile 
import pandas as pd
import os
import csv


# In[2]:


def openzip(ZIPFILE, time_resolution):
    
    '''Input: ZIPFILE, time_resolution(YY: Yearly, MM: Monthly, DD: Daily, WW: weekly, HH: hourly);
    Output: df: dataframe, name: name of data file;
    Opens ZIPFILE and searches for given time resolution and returns the name of the file and the datafile'''

    zf = zipfile.ZipFile(ZIPFILE)
    names = zf.namelist()
    for i in range(len(names)):
        if (names[i].find('FULLSET_' + time_resolution) != -1):
            print (str(names[i])+ " found")
            name = names[i]
            df = pd.read_csv(zf.open(name))
            
    return df,name

def plot_data(file, time_resolution,Index,path,path_goal):
    
    ''' Input: file, time_resolution(YY: Yearly, MM: Monthly, DD: Daily, WW: weekly, HH: hourly), Index: wanted parameter, 
    path: path of data file, path_goal: path to save figure;
    Selects variables in data frame with given time_resolution, plots figure in path_goal
    '''
    
    ZIPFILE = path + file
    df, name = openzip(ZIPFILE, time_resolution)
    
    for i in range(len(Index)):
        
        path_ = path_goal + Index[i]
        #select wanted variable
        data_ = df[["TIMESTAMP",Index[i]]]
        #filter NaN-values
        data__ = data_[data_[Index[i]] > - 9999]
        
        plt.figure()
        plt.plot(data__["TIMESTAMP"],data__[Index[i]],'*')
        plt.title(name[:-3])
        plt.xlabel(time_resolution)
        plt.ylabel(Index[i])
        
        #create folder for variable
        try: 
            os.mkdir(path_)
        except:
            a = 0
            
        plt.savefig(path_ + "/" + str(name[:-3]) + "_" + time_resolution + ".png")
        
        
def select_data(file,time_resolution, Index, Index_Names, path,path_goal):
    '''Input: file, time_resolution(YY: Yearly, MM: Monthly, DD: Daily, WW: weekly, HH: hourly), Index: wanted parameters, 
    Index_names: names for new dataframe, path: path of data file, path_goal: path to save csv file;
    searches for wanted parameters and creates a new data frame with reduced parameter. The new file is than stored in path_goal'''
    
    ZIPFILE = path + file
    df,name = openzip(ZIPFILE, time_resolution)
    data = df[Index]
    data.columns = Index_Names
    file_goal = path_goal+file[:-4]+"_reduced.csv"
    data.to_csv(file_goal)
    


# In[ ]:




