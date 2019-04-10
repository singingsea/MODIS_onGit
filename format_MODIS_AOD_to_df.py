# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 11:16:02 2018

@author: xiaoy
"""


#!/usr/bin/python
'''
Module: read_mod_aerosol_at_a_location.py
==========================================================================================
Disclaimer: The code is for demonstration purposes only. Users are responsible to check for accuracy and revise to fit their objective.

Author: Justin Roberts-Pierel, 2015 
Organization: NASA ARSET
Purpose: To view info about a variety of SDS from a MODIS HDF4 file (or series of files) both generally and at a specific lat/lon

See the README associated with this module for more information.
==========================================================================================
'''

#import necessary modules
from pyhdf import SD
import numpy as np
import pandas as pd
import os


#MODIS_data_path = 'E:\\Projects\\MODIS\\501219332\\'
MODIS_data_path = 'C:\\Projects\\MODIS\\data\\'
#output_file = 'E:\\Projects\\MODIS\\test.csv'
output_file = ':\\Projects\\MODIS\\output\\test.csv'
user_lat = 40
user_lon = -50
#user_lat = 43
#user_lon = -73
os.chdir(MODIS_data_path)
df = pd.DataFrame()

#This uses the file "fileList.txt", containing the list of files, in order to read the files
try:
	fileList = open(str(MODIS_data_path + 'fileList.txt'),'r')
except:
	print('Did not find a text file containing file names (perhaps name does not match)')
	sys.exit()

dataFields=dict([(1,'Optical_Depth_Land_And_Ocean'),(2,'Land_Ocean_Quality_Flag'), (3,'Image_Optical_Depth_Land_And_Ocean'),(4,'Land_Sea_Flag')])
userInput = 1    
#loops through all files listed in the text file
i=0
for FILE_NAME in fileList:
    FILE_NAME=FILE_NAME.strip()

    if '3K' in FILE_NAME: #then this is a 3km MODIS file
        dataFields =dict([(1,'Optical_Depth_Land_And_Ocean'),(2,'Land_Ocean_Quality_Flag'),(3,'Image_Optical_Depth_Land_And_Ocean'),(4,'Land_Sea_Flag')])
    elif 'L2' in FILE_NAME:#Same as above but for 10km MODIS file
        dataFields=dict([(1,'Deep_Blue_Aerosol_Optical_Depth_550_Land'),(2,'AOD_550_Dark_Target_Deep_Blue_Combined'),(3,'AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag')])
	
    SDS_NAME=dataFields[int(userInput)] # The name of the sds to read

    # open the hdf file for reading
    hdf=SD.SD(FILE_NAME)

		
    # Get lat and lon info
    lat = hdf.select('Latitude')
    latitude = lat[:,:]
    min_lat=latitude.min()
    max_lat=latitude.max()
    lon = hdf.select('Longitude')
    longitude = lon[:,:]
    min_lon=longitude.min()
    max_lon=longitude.max()
		
    #get SDS, or exit program if SDS is not in the file
    try:
        sds=hdf.select(SDS_NAME)
    except:
        print('Sorry, your MODIS hdf file does not contain the SDS:',SDS_NAME,'. Please try again with the correct file type.')
        continue
    #get scale factor and fill value for data field
    attributes=sds.attributes()
    scale_factor=attributes['scale_factor']
    fillvalue=attributes['_FillValue']
    
    #get SDS data
    data=sds.get()
    #Print the range of latitude and longitude found in the file, then ask for a lat and lon
    print('The range of latitude in this file is: ',min_lat,' to ',max_lat, 'degrees \nThe range of longitude in this file is: ',min_lon, ' to ',max_lon,' degrees')
    #user_lat=float(input('\nPlease enter the latitude you would like to analyze (Deg. N): '))
    #user_lon=float(input('Please enter the longitude you would like to analyze (Deg. E): '))
    #Continues to ask for lat and lon until the user enters valid values
    escape = False
    while user_lat < min_lat or user_lat > max_lat:
        print('The latitude you entered is out of range. Please check your data file! ')
        escape = True
    while user_lon < min_lon or user_lon > max_lon:
        print('The longitude you entered is out of range. Please check your data file! ')
        escape = True
    if escape != True:
        #calculation to find nearest point in data to entered location (haversine formula)
        R=6371000#radius of the earth in meters
        lat1=np.radians(user_lat)
        lat2=np.radians(latitude)
        delta_lat=np.radians(latitude-user_lat)
        delta_lon=np.radians(longitude-user_lon)
        a=(np.sin(delta_lat/2))*(np.sin(delta_lat/2))+(np.cos(lat1))*(np.cos(lat2))*(np.sin(delta_lon/2))*(np.sin(delta_lon/2))
        c=2*np.arctan2(np.sqrt(a),np.sqrt(1-a))
        d=R*c
        #gets (and then prints) the x,y location of the nearest point in data to entered location, accounting for no data values
        x,y=np.unravel_index(d.argmin(),d.shape)
        print('\nThe nearest pixel to your entered location is at: \nLatitude:',latitude[x,y],' Longitude:',longitude[x,y])
        if data[x,y]==fillvalue:
            print('The value of ',SDS_NAME,'at this pixel is',fillvalue,',(No Value)\n')
        else:
            print('The value of ', SDS_NAME,'at this pixel is ',round(data[x,y]*scale_factor,3))
        
        var_nm = str(SDS_NAME +'_nearest_pixel')
        kw = {var_nm:'NaN'}
        df = df.assign(**kw)
        #df[var_nm] = 0
        df[var_nm].iloc[0:i]	= round(data[x,y]*scale_factor,3)
        i += 1
        #calculates mean, median, stdev in a 3x3 grid around nearest point to entered location
        if x < 1:
            x+=1
        if x > data.shape[0]-2:
            x-=2
        if y < 1:
            y+=1
        if y > data.shape[1]-2:
            y-=2
        three_by_three=data[x-1:x+2,y-1:y+2]
        three_by_three=three_by_three.astype(float)
        three_by_three[three_by_three==float(fillvalue)]=np.nan
        nnan=np.count_nonzero(~np.isnan(three_by_three))
        if nnan == 0:
            print ('\nThere are no valid pixels in a 3x3 grid centered at your entered location.')
        else:
            three_by_three=three_by_three*scale_factor
            three_by_three_average=np.nanmean(three_by_three)
            three_by_three_std=np.nanstd(three_by_three)
            three_by_three_median=np.nanmedian(three_by_three)
            if nnan == 1:
                npixels='is'
                mpixels='pixel'
            else:
                npixels='are'
                mpixels='pixels'
            print('\nThere',npixels,nnan,'valid',mpixels,'in a 3x3 grid centered at your entered location.')
            print('\nThe average value in this grid is: ',round(three_by_three_average,3),' \nThe median value in this grid is: ',round(three_by_three_median,3),'\nThe standard deviation in this grid is: ',round(three_by_three_std,3))
        
    #    df[str(SDS_NAME +'_3x3_pixel')].loc(i)	= round(three_by_three_average,3)
        
        #calculates mean, median, stdev in a 5x5 grid around nearest point to entered location
        if x < 2:
            x+=1
        if x > data.shape[0]-3:
            x-=1
        if y < 2:
            y+=1
        if y > data.shape[1]-3:
            y-=1
            five_by_five=data[x-2:x+3,y-2:y+3]
            five_by_five=five_by_five.astype(float)
            five_by_five[five_by_five==float(fillvalue)]=np.nan
            nnan=np.count_nonzero(~np.isnan(five_by_five))
        if nnan == 0:
            print ('There are no valid pixels in a 5x5 grid centered at your entered location. \n')
        else:
            five_by_five=five_by_five*scale_factor
            five_by_five_average=np.nanmean(five_by_five)
            five_by_five_std=np.nanstd(five_by_five)
            five_by_five_median=np.nanmedian(five_by_five)
            if nnan == 1:
                npixels='is'
                mpixels='pixel'
            else:
                npixels='are'
                mpixels='pixels'
            print('\nThere',npixels,nnan,' valid',mpixels,' in a 5x5 grid centered at your entered location. \n')
            print('The average value in this grid is: ',round(five_by_five_average,3),' \nThe median value in this grid is: ',round(five_by_five_median,3),'\nThe standard deviation in this grid is: ',round(five_by_five_std,3))       
        var_nm = str(SDS_NAME +'_5x5_pixel')
#    df[var_nm].loc(i) = round(five_by_five_average,3)  
        
print('\nAll valid files have been processed')
df.to_csv(output_file)