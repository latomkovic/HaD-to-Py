from hec.script import *
from hec.heclib.dss import *
from hec.hecmath import *
from hec.heclib.util import HecTime
import java
import sys, hec, time, os
import os.path
from hec.io import TimeSeriesContainer

# Since spaces separate each index of the sys.argv list, this puts together a path string
def get_dss_filename(i):
	dss_ext = False
	filename_list = []
	if arg_list[i][-3:].lower() == 'dss':
		filename = arg_list[i]
		filename = ''.join(filename)
		dss_ext = True
		i+=1
	while dss_ext == False:
		filename_list.append(arg_list[i])
		filename = ' '.join(filename_list)
		if arg_list[i][-3:].lower() == 'dss':
			dss_ext = True
		i+=1
	return filename, i

def get_txt_filename(i):
	txt_ext = False
	filename_list = []
	if arg_list[i][-3:].lower() == 'txt':
		filename = arg_list[i]
		filename = ''.join(filename)
		txt_ext = True
		i+=1
	while txt_ext == False:
		filename_list.append(arg_list[i])
		filename = ' '.join(filename_list)
		if arg_list[i][-3:].lower() == 'txt':
			txt_ext = True
		i+=1
	return filename	

def get_path_from_sys(i):
	slash = 0
	path_list = []
	while (slash < 7):
		path_list.append(arg_list[i])
		path = ''.join(path_list)
		if path[-1:] != '/':
			path_list.append(' ')
		slash = path.count("/")
		i+=1
	return path, i

def get_datetime_from_sys(i):
	time = []
	time.append(arg_list[i])
	time.append(' ')
	time.append(arg_list[i+1])
	datetime = ''.join(time)
	i +=2
	return datetime,i
	
def get_data_list(i):
	data_count = int(arg_list[i])
	i+=1
	count = 0
	data_list = []
	while count < data_count:
		data_list.append(float(arg_list[i]))
		count+=1
		i+=1
	return data_list
	
# Initialize argument list from system command input
arg_list = []
for arg in sys.argv:
	arg_list.append(arg) 

dss_filename, data_ind = get_dss_filename(1)
#MessageBox.showPlain("filename: %s"%dss_filename, "test")
data_units = arg_list[data_ind]
#MessageBox.showPlain("DataUnits: %s"%data_units, "test")
pathname, start_ind = get_path_from_sys(data_ind+1)
#MessageBox.showPlain("pathname: %s"%pathname, "test")
start_time, comp_ind =get_datetime_from_sys(start_ind)
#MessageBox.showPlain("start_time: %s"%start_time, "test")
comp_step = int(arg_list[comp_ind])
#MessageBox.showPlain("computation step: %s"%comp_step, "test")
data_count = int(arg_list[comp_ind+1])
input_file = get_txt_filename(comp_ind+2)
#MessageBox.showPlain("data list: %s"%input_file, "test")

if not os.path.isfile(dss_filename):
	MessageBox.showError("DSS was unable to open the DSS file. \n The filename provided is: %s \n Please check to see that your RAS project is exporting data to this file!"%(dss_filename),"DSS unable to open file!")

# Read from input_file and produce a list to store to TS record
data_list = []	
txt_input = open(input_file,'r')
for line in txt_input:
	data_list.append(float(line))
txt_input.close()


def createTSrecord(dss_filename, pathname, start_time, values, comp_step, data_units ):
	
	
	start = HecTime()
	tsc = TimeSeriesContainer()
	tsc.fullName = pathname
	tsc.interval = comp_step
	start.set(start_time)
	times =  []
	for value in values:
		times.append(start.value())
		start.add(tsc.interval)
	tsc.values = values
	tsc.times = times
	tsc.startTime = times[0]
	tsc.endTime = times[-1]
	tsc.numberValues = len(values)
	tsc.units = data_units	
	tsc.type = "INST-VAL"
	dss_file = HecDss.open(dss_filename)
	dss_file.put(tsc) 
	dss_file.done()



createTSrecord(dss_filename, pathname, start_time, data_list, comp_step, data_units )










