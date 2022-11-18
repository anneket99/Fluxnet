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
    data.to_csv(file_goal, index = False)
    

def location(file,path,path_goal):
    
    ''' Input: file ( with information on location), path (with fluxnet files), path_goal (path for storing data)
    
    
    The code first reads in the location and SITE_ID of the input file and then looks for the fluxnet sites with the same SITE_ID. 
    The SITE_ID and loccation get added to the dataset. All datasets in the given dataset get combined into one dataset.'''
    
    result = []
    
    # read file with locations
    with open(file) as f:
        f = open(file,"r")
        for x in f:
            result.append(x.split(';')[1])
            
    SITE_ID_ = [x.split(';')[0] for x in open(file).readlines()]
    SITE_ID = []
    for i in range (1,len(SITE_ID_)):
        SITE_ID.append(SITE_ID_[i][1:])
    
    #get latitute and longitude
    LOCATION_LAT = [x.split(';')[4] for x in open(file).readlines()][1:]
    LOCATION_LONG = [x.split(';')[5] for x in open(file).readlines()][1:]
    
    #look for files with same SITE_ID 
    names = os.listdir(path)
    data = []
    for j in range(len(names)):
        if (names[j].find('reduced') != -1):
            file = path+names[j]
            data_ = pd.read_csv(file)
        
            for i in range(len(SITE_ID)):
                
                # add SITE_ID, Latitude and Longitude
                if(file.find(SITE_ID[i]) != -1):
                    data_.insert(1,"SITE_ID",SITE_ID[i])
                    data_.insert(2,"LAT",LOCATION_LAT[i])
                    data_.insert(3,"LONG", LOCATION_LONG[i])
                    try:
                        data = data.append(data_,ignore_index = True)
                    except:
                        data = data_
    data.to_csv(path_goal, Index = False)




