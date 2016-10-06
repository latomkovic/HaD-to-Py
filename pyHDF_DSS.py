#!/usr/bin/env python
"""
Takes HDF and DSS files (outputs from HEC-RAS) and puts together plots of 
observed versus computed water surface elevations at user-specified gage
locations.
 
Created on Mon Aug 22 16:25:09 2016

@author: lelekew

"""

import os
import numpy as np
from scipy import spatial
from Tkinter import *
import tkMessageBox
import subprocess
import h5py 
import matplotlib.pyplot as plt

# User input:
hdf_filename = r"C:/Users/lelekew/Documents/Python/Example RAS/Example Projects/2D Unsteady Flow Hydraulics/BaldEagleCrkMulti2D/BaldEagleDamBrk.p03.hdf" 
obs_dss = r'C:/Users/lelekew/Documents/Python/FinalFolder/Tutorial/BaldCreekObservedData.dss' 
twoD_dss = r"C:/Users/lelekew/Documents/Python/FinalFolder/Tutorial/BaldCreek2D.dss" 
plot_dir = r'C:\Users\lelekew\Documents\Python\Figures' 

working_directory = os.path.dirname(os.path.realpath(__file__))

def make_DSS_plot(plan, gage, ras_dss, obs_dss, plot_dir, ras_path, obs_path, start_time, \
    end_time, run_time):
    path = r"C:\\Program Files (x86)\\HEC\\HEC-DSSVue"
    executable = "HEC-DSSVue.exe "
    if not os.path.exists(path):
        message = """Script needs HEC-DSSVue in order to run. 
            Install HEC-DSSVue and place in C:\Program Files (x86)\HEC\HEC-DSSVue
            or Change path variable in make_DSS_plot function.
            """
        error_box(message)
    
    infile = working_directory + "\\DSSplot.py "
    if not os.path.isfile(infile):
        message = """Script needs getDSSdata.py file to be in the working directory.
            The specified file is in the download package for this script.
            """
        error_box(message)
    
    cmds = ['pushd "%s" &&' % path, executable, infile, plan, gage, ras_dss, obs_dss, \
            plot_dir, ras_path, obs_path, start_time, end_time, run_time, "&& popd"]
    
    subprocess.call(
                " ".join(cmds),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True)

# Error boxes for future
def error_box(message_text):
    window = Tk()
    window.wm_withdraw()
    window.geometry("1x1+"+str(window.winfo_screenwidth()/2)+"+"+str(window.winfo_screenheight()/2))
    tkMessageBox.showinfo(title="Error", message=message_text)
    window.lift()
    window.mainloop()
    

# Get the short identifier    

def get_shortID(hdf_filename):
    plan_filename = hdf_filename[:-4]
    ras_dss = plan_filename[:-4] + '.dss'
    plan = plan_filename[-3:]
    short_indicator = "Short Identifier="
    project_filename = plan_filename[:-4] + '.prj'
    with open(project_filename,'r') as project_file:
        lines = project_file.readlines()
        if "English Units" in lines[3]:
            units = 'Feet'
        elif "SI Units" in lines[3]:
            units = "Meter"
        else:
            units = "Somehow Unknown"
    with open(plan_filename,'r') as plan_file:
        for line in plan_file:
            if short_indicator in line:
                shortID = line[line.index(short_indicator) + len(short_indicator):].rpartition("\n")[0]
                shortID = ' '.join(shortID.split())
                break
    return shortID, plan, ras_dss, units  
    
# Get the computation intervals
def get_intervals(hdf_filename):
    plan_filename = hdf_filename[:-4]
    comp_indicator = "Computation Interval="
    outp_indicator = "Output Interval="
    map_indicator = "Mapping Interval="
    i=0
    with open(plan_filename,'r') as plan_file:
        for line in plan_file:
            if comp_indicator in line:
                comp_interval = line[line.index(comp_indicator) + len(comp_indicator):].rpartition("\n")[0]
                i += 1
                if i==3:
                    break
            elif outp_indicator in line:
                outp_interval = line[line.index(outp_indicator) + len(outp_indicator):].rpartition("\n")[0]
                i += 1
                if i==3:
                    break
            elif map_indicator in line:
                map_interval = line[line.index(map_indicator) + len(map_indicator):].rpartition("\n")[0]
                i += 1
                if i==3:
                    break
    return comp_interval, outp_interval, map_interval    

# Pass the interval to this function to get the minute value
def get_time_minutes(interval):
    interval_to_minutes = {'1MIN': 1, '2MIN':2,'3MIN':3,'4MIN':4,'5MIN':5,'10MIN':10,'15MIN':15,'20MIN':20,'30MIN':30,'1HOUR':60}
    try:
        minutes = interval_to_minutes[interval]
    except:
        message = "You might need to modify the get_time_minutes() function so that there is a interval corresponding to your %s interval" %interval
        error_box(message)
        raise
    return minutes   
   
# Get the 2D Coordinates from the txt file
def get_coordinates(working_directory):
    coords_txt_file = working_directory+"\\two_dim_coords.txt"
    coordinates = {}
    if os.path.isfile(coords_txt_file):
        with open(coords_txt_file,'r') as coords_file:
            for line in coords_file:
                key = line.partition('Gage: ')[-1].rpartition(' Lat:')[0]
                latitude = line.partition('Lat: ')[-1].rpartition(' Lon')[0]
                longitude = line.partition('Lon: ')[-1]
                latitude = float(latitude)
                longitude = float(longitude)
                coordinates.setdefault(key,[]).append(latitude)
                coordinates.setdefault(key,[]).append(longitude)
                
    else:        
        message="""The file 'two_dim_coords.txt' needs to be filled with coordinates and placed in the working directory.
                The format is:
                Gage: Gagename Lat: 123456.789 Lon: 987654.321
                ex:
                Gage: FRE Lat: 615795.406826 Lon: 4290904.93965
                Where the coordinates are latitude and longitude in meters."""
        error_box(message)
    
    return coordinates
        

# Get the plan_filename and the start_time and end_time
def get_start_end_time(hdf_filename):
    plan_filename = hdf_filename[:-4]
    time_indicator = "Simulation Date="
    with open(plan_filename,'r') as plan_file:
        for line in plan_file:
            if time_indicator in line:
                time_line = line
                break
    times = time_line[time_line.index(time_indicator) + len(time_indicator):]
    start_time = times[0:14]
    start_time = start_time.replace(","," ")
    end_time = times[15:]
    end_time = end_time.replace(","," ")
    Dpart = start_time[0:9]
    Dpart = list(Dpart)
    Dpart[0:2] = "01"
    Dpart = ''.join(Dpart)
    return start_time, end_time, Dpart   
    
# Get the observed data path from the gage name
def get_obs_path(gage):
    obs_pathname_file = working_directory + "\\obs_paths.txt"
    with open(obs_pathname_file,'r') as obs_file:
        for line in obs_file:
            if gage in line:
                obs_path = line.partition('Path: ')[-1].rpartition('\n')[0]
    try:
        obs_path
    except UnboundLocalError:
        error_box("Please ensure that there is a DSS path corresponding to the observed DSS file in the 'one_dim_comp_paths' found in the working directory.")
        raise
    return obs_path
    
# Simply gets the names of the 2D Flow Areas in the Plan's geometry    
def get2DAreaNames(hf):
    hdf2DFlow = hf['Geometry']['2D Flow Areas']
    AreaNames = []
    for key in hdf2DFlow:
        if key in ['Cell Spacing', "Manning's n", 'Names', 'Tolerances']:
            continue
        else:
            AreaNames.append(key) # List of 2D Area names
    return AreaNames

# Gets the nearest cell in a given area, needs to be run through findGageCell in order to find absolute closest cell. 
# (i.e. it could return the closest cell in Area A but the gage is actually in Area B)
def findNearestCell(hf,AreaName,gageCoord):
    # Finds the nearest cell to the gageCoordinate from the HDF file
    
    CellPts = hf['Geometry']['2D Flow Areas'][AreaName]['Cells Center Coordinate'][:] # All Coordinates which are in 2D Flow Area with name AreaName
    fPtInd = hf['Geometry']['2D Flow Areas'][AreaName]['Cells FacePoint Indexes'] # All FacePoint Indices which correspond to each 'cell center point'
    # the reason this next step is necessary is that CellPts returns coordinates for face points as well as center points, and we need the index of the cell point.
    distances,indices = spatial.KDTree(CellPts).query(gageCoord, k=7) # finds the closest points (cell center or face)
    # Now to eliminate FacePoints and only have cell centers left
    i = 0
    for cell in range(0, len(indices)):
        if fPtInd[indices[i]][2] == -1:
            distances=np.delete(distances,i)
            indices=np.delete(indices,i)
        else:
            i+=1
            continue
    if indices.size: # If the indices list is not empty
        CellInd = indices[np.where(distances==min(distances))] # Choose the index that's left with the shortest distance
        CellDist = distances[np.where(distances==min(distances))]  # Displays the distance left
    else:
        CellInd = [0] # Be careful of this
        CellDist = [9999] # Seeing this value in the Distances array will notify you that none of the cells in the 2D flow area were associated with a cell point. 
                        # This is likely because the gage coordinates are too far from a 2D area to be considered "close" to a cell point 
                        # and so face points are all that are being rendered as "close"
    return CellInd[0], CellDist[0]  # The index of the cell center which is closest to the gage coordinates

# Finds the absolute closest cell to the gage coordinates given. Returns the cell index to be used to access
# output data and the distance of that cell's center to the gage coordinates.
def findGageCell(hf,coord):
    AreaNames = get2DAreaNames(hf)
    cell=list()
    dist=list()
    areas=list()
    for area in AreaNames:
        tempInd,tempDist = findNearestCell(hf,area,coord)
        cell.append(tempInd)
        dist.append(tempDist)
        areas.append(area)
    cellInd = cell[dist.index(min(dist))] # Choose the index that's left with the shortest distance
    cellDist = dist[dist.index(min(dist))]  # Displays the distance left
    cellArea = areas[dist.index(min(dist))]
    return cellInd, cellDist, cellArea

# Uses the coordinates in the 'two_dim_coords.txt' file in the working directory to access the cell indices of the
# 2D areas in the hdf_file specified by the variable hdf_filename. Returns the dictionaries Distances and Indices
# as well as a txt file which contains that information so as to be accessed later in a temp file
def find_2D_indices(working_directory,hdf_filename):
    coordinates = get_coordinates(working_directory)
    Indices = coordinates.fromkeys(coordinates,0)
    Distances = coordinates.fromkeys(coordinates,0)
    Areas = coordinates.fromkeys(coordinates,0)
    hf = h5py.File(hdf_filename,'r') 
    
    # Iterates through all of the gages given in the coordinates dictionary to find and store the distances 
    # and cell indices corresponding to each gage.
    indices_outfile = working_directory + "\\temp_files\\ras_2D_cells.txt"
    text_file = open(indices_outfile, "w")
    text_file.write(hdf_filename + '\n')
    for gage in coordinates:
        cellInd, cellDist, cellArea = findGageCell(hf,coordinates[gage])
        Distances.update({gage:cellDist})
        Indices.update({gage:cellInd})
        Areas.update({gage:cellArea})
        text_file.write("Gage: %s Distance: %s Area: %s Index: %s\n" %(gage,Distances[gage],Areas[gage],Indices[gage]))
    text_file.close() 
    return Distances, Indices, Areas

# Check to see if a 
def get_2D_indices(working_directory,hdf_filename):
    temp_file_dir = working_directory + "\\temp_files"
    if not os.path.exists(temp_file_dir): os.makedirs(temp_file_dir)
    if not os.path.isfile(temp_file_dir + "\\ras_2D_cells.txt"): indices_file= open(temp_file_dir+"\\ras_2D_cells.txt",'w')
    indices_file = open(working_directory + "\\temp_files\\ras_2D_cells.txt", 'r+')
    lines = indices_file.readlines()
    indices_file.close()
    if len(lines)<2:
        headline = ''
    else:
        headline = lines[0]
        iter_lines = lines[1:]
    if hdf_filename in headline:
        #If the ras_2D_cells.txt file is already populated with the current plan,
        #then the script doens't need to re-find the indices etc.
        Distances = {}
        Indices = {}
        Areas = {}
        for line in iter_lines:
            key = line.partition('Gage: ')[-1].rpartition(' Distance:')[0]
            distance = line.partition('Distance: ')[-1].rpartition(' Area:')[0]
            area = line.partition('Area: ')[-1].rpartition(' Index:')[0]
            index = line.partition('Index: ')[-1].rpartition('\n')[0]
            Distances[key]=float(distance)
            Indices[key]=int(index)
            Areas[key]=str(area)
    else:
        #Find the cell indices for the coordinates given in the two_dim_coords.txt file
        Distances, Indices, Areas = find_2D_indices(working_directory,hdf_filename)
    return Distances, Indices, Areas

def get_DSS_data(dss_filename, temp_output_filename, pathname,start_time, end_time):
    path = r"C:\\Program Files (x86)\\HEC\\HEC-DSSVue"
    executable = "HEC-DSSVue.exe "
    
    if not os.path.exists(path):
        message = """Script needs HEC-DSSVue in order to run. 
            Install HEC-DSSVue and place in C:\Program Files (x86)\HEC\HEC-DSSVue
            or Change path variable in make_DSS_plot function.
            """
        error_box(message)
    
    infile = working_directory + "\\getDSSdata.py "
    
    if not os.path.isfile(infile):
        message = """Script needs getDSSdata.py file to be in the working directory.
            The specified file is in the download package for this script.
            """
        error_box(message)
    
    cmds = ['pushd "%s" &&' % path, executable, infile, dss_filename, \
            temp_output_filename, pathname, start_time, end_time,  "&& popd"]   
    
    subprocess.call(
                " ".join(cmds),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True)
    
    # HEC-DSSVue has created a temp file (temp_output_filename) and now we read it into a data list
    temp_data_file = open(temp_output_filename, 'r+')
    lines = temp_data_file.readlines()
    if "NO DATA" in lines[0]:
        data_list = "NO DATA"
        time_str_list = "NO DATA"
        time_min_list = "NO DATA"
    else:
        data_line = lines[1]
        data_line = data_line.replace("array('d',[","")
        data_line = data_line.replace('])','')
        data_list = data_line.split(",")
        data_list = [float(i) for i in data_list]
        time_str_line = lines[3]
        time_str_line = time_str_line.replace("[","")
        time_str_line = time_str_line.replace(']','')
        time_str_list = time_str_line.split("',")
        for t in range(len(time_str_list)):
            time_str_list[t] = time_str_list[t].replace("'",'')
        time_min_line = lines[5]
        time_min_line = time_min_line.replace("array('i',[","")
        time_min_line = time_min_line.replace('])','')
        time_min_list = time_min_line.split(",")
        time_min_list = [int(i) for i in time_min_list]
    temp_data_file.close()
    return data_list, time_str_list, time_min_list
    
# Store to DSS
def store_to_DSS(data_units, temp_input_filename, dss_filename, pathname, start_time, comp_step, data_list):
    # comp_step is an integer value representing the number of minutes in a timestep (1HOUR - 60, etc)
    
    path = r"C:\\Program Files (x86)\\HEC\\HEC-DSSVue"
    executable = "HEC-DSSVue.exe "
    
    if not os.path.exists(path):
        message = """Script needs HEC-DSSVue in order to run. 
            Install HEC-DSSVue and place in C:\Program Files (x86)\HEC\HEC-DSSVue
            or Change path variable in make_DSS_plot function.
            """
        error_box(message)
    
    infile = working_directory + "\\storeDSSdata.py "
    comp_step = str(comp_step)
    data_units = data_units.replace(' ', '')   
    data_count = str(len(data_list))
    #data_list = ' '.join(str(x) for x in data_list)
    
    # HEC-DSSVue maxes out at 100 command-line entries, which is very limiting for data entry! So must read from txt file.
    with open(temp_input_filename, 'w') as input_file:
        for entry in data_list:        
            input_file.write("%s \n" %str(entry))
        input_file.close()
        
    if not os.path.isfile(infile):
        message = """Script needs getDSSdata.py file to be in the working directory.
            The specified file is in the download package for this script.
            """
        error_box(message)
    
    cmds = ['pushd "%s" &&' % path, executable, infile, dss_filename, data_units,\
            pathname, start_time, comp_step, data_count, temp_input_filename,  "&& popd"]
    
    subprocess.call(
                ' '.join(cmds),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True)
                
def get_HDF_cell_WSE(hf, cell_number, flow_area):
    
    with h5py.File(hdf_filename,'r') as hf:
        
        flow_areas = hf['Results']['Unsteady']['Output']['Output Blocks']\
        ['Base Output']['Unsteady Time Series']['2D Flow Areas']
        
        dataset = flow_areas[flow_area]['Water Surface']
        timesteps = dataset.shape[0]
        
        data_list = np.zeros((timesteps,), dtype='Float64')
        dataset.read_direct(data_list, np.s_[0:timesteps,cell_number], np.s_[0:timesteps])
        data_list = np.array(data_list).tolist()
            
    return data_list                            


# This will go through all of the 1D and 2D observed points listed in the two_dim_coords and one_dim_comp_paths txt files
# Without those two files, the program will not run. This function returns data dictionaries for each gage            
def get_all_plot_data(working_directory, obs_dss, twoD_dss, hdf_filename):
    Distances, Indices, Areas = get_2D_indices(working_directory,hdf_filename)
    shortID, plan, ras_dss, units = get_shortID(hdf_filename)
    temp_output_filename = working_directory + "\\temp_files\\DSS_data_output.txt"
    temp_input_filename = working_directory + "\\temp_files\\DSS_data_input.txt"
    start_time, end_time, Dpart = get_start_end_time(hdf_filename)
    comp_interval, outp_interval, map_interval = get_intervals(hdf_filename)
    map_minutes = get_time_minutes(map_interval)
    
    # Initiate dictionaries which will contain the data for each gage location
    ras_computed_WSE = {}
    observed_WSE = {}
    time_strings_obs = {}
    time_strings_comp = {}
    time_minutes_obs = {}
    time_minutes_comp = {}
    
    # Get all 2D Data
    for key in Indices:
        pathname = "/%s//STAGE/%s/%s/%s/" %(key,Dpart,map_interval,shortID)
        ras_computed_WSE[key] = get_HDF_cell_WSE(hdf_filename, Indices[key], Areas[key])
        store_to_DSS(units, temp_input_filename, twoD_dss, pathname, start_time, map_minutes, ras_computed_WSE[key])
        obs_path = get_obs_path(key)
        observed_WSE[key], time_strings_obs[key], time_minutes_obs[key] = get_DSS_data(obs_dss, temp_output_filename, obs_path,start_time, end_time)
        time_minutes_comp[key]=range(len(ras_computed_WSE[key]))
        map_minutes = get_time_minutes(map_interval)
        for time in range(len(ras_computed_WSE[key])):        
            if time == 0:
                time_minutes_comp[key][0] = time_minutes_obs[key][0] 
            else:
                time_minutes_comp[key][time] = time_minutes_comp[key][time-1] + map_minutes
        
            
    # Get all the 1D Data
    if os.path.isfile(working_directory + "\one_dim_comp_paths.txt"):
        one_dim_comp_path_file = open(working_directory + "\one_dim_comp_paths.txt", 'r+')
        for line in one_dim_comp_path_file:
            key = line.partition('Gage: ')[-1].rpartition(' Path:')[0]
            comp_path = line.partition('Path: ')[-1].rpartition('\n')[0]
            comp_path = comp_path %(Dpart, outp_interval, shortID)
            obs_path = get_obs_path(key)
            ras_computed_WSE[key], time_strings_comp[key], time_minutes_comp[key] = get_DSS_data(ras_dss, temp_output_filename, comp_path, start_time, end_time)
            observed_WSE[key], time_strings_obs[key], time_minutes_obs[key] = get_DSS_data(obs_dss, temp_output_filename, obs_path,start_time, end_time)
    one_dim_comp_path_file.close()
    
    return ras_computed_WSE, observed_WSE, time_strings_comp, time_minutes_comp, time_strings_obs, time_minutes_obs

# Plot observed and computed and save to plot directory
def plot_and_save(time1,value1, time2,value2, plot_name, plot_dir):
    if time1 == "NO DATA":
        plt.plot(time2,value2,'r-', label='Computed')
        plt.ylabel('Stage (m)')
        plt.xlabel('Minutes HEC Time')
        plt.legend( loc = 'upper left', numpoints = 1)
        plt.title(plot_name)
        plt.savefig(plot_dir + '\\' + plot_name + '.png')
        plt.close()
    elif time2 == "NO DATA":
        plt.plot(time1,value1,'b', linestyle = '-',marker = 'o', markersize = 3, label='Observed')
        plt.ylabel('Stage (m)')
        plt.xlabel('Minutes HEC Time')
        plt.legend( loc = 'upper left', numpoints = 1)
        plt.title(plot_name)
        plt.savefig(plot_dir + '\\' + plot_name + '.png')
        plt.close()
    else:
        plt.plot(time1,value1,'b', linestyle = '-',marker = 'o', markersize = 3, label='Observed')
        plt.plot(time2,value2,'r-', label='Computed')
        plt.ylabel('Stage (m)')
        plt.xlabel('Minutes HEC Time')
        plt.legend( loc = 'upper left', numpoints = 1)
        plt.title(plot_name)
        plt.savefig(plot_dir + '\\' + plot_name + '.png')
        plt.close()

# Go through all gages and plot them and save them in a folder under the shortID plan name
def store_plots(obs_time, obs_value, comp_time, comp_value, plot_directory, hdf_filename):
    shortID, plan, ras_dss, units = get_shortID(hdf_filename)
    plot_subfolder = plot_directory + '\\' + shortID
    if not os.path.exists(plot_subfolder): 
        os.makedirs(plot_subfolder)
    for key in comp_value:
        plot_and_save(obs_time[key],obs_value[key], comp_time[key],comp_value[key], key, plot_subfolder)
        
ras_computed_WSE, observed_WSE, time_strings_comp, time_minutes_comp, time_strings_obs, time_minutes_obs = get_all_plot_data(working_directory, obs_dss, twoD_dss, hdf_filename)
store_plots(time_minutes_obs, observed_WSE, time_minutes_comp, ras_computed_WSE, plot_dir, hdf_filename)

'''
# Here lie some function calls
outputlist = get_HDF_cell_WSE(hdf_filename, Indices['FRE'], Areas['FRE'])
Distances, Indices, Areas = get_2D_indices(working_directory,hdf_filename)
store_to_DSS(data_units, test_dss, test_path, start_time, comp_step, test_data_list)
shortID, plan, ras_dss, units = get_shortID(hdf_filename)
comp_interval, outp_interval, map_interval = get_intervals(hdf_filename)
start_time, end_time, Dpart = get_start_end_time(hdf_filename)
data_list, time_str_list, time_min_list = get_DSS_data(obs_dss, temp_output_filename, obs_path,start_time, end_time)
'''
