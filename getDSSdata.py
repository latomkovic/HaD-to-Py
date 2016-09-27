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
	return filename, i
		

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

# Initialize argument list from system command input
arg_list = []
for arg in sys.argv:
	arg_list.append(arg) 
#MessageBox.showPlain("%s" %arg_list,"Test")

dss_filename, output_ind = get_dss_filename(1)
#MessageBox.showPlain("%s" %dss_filename,"Test")
output_filename, path_ind = get_txt_filename(output_ind)
#MessageBox.showPlain("%s" %output_filename,"Test")
pathname, start_ind = get_path_from_sys(path_ind)
#MessageBox.showPlain("%s" %pathname,"Test")
start_time, end_ind =get_datetime_from_sys(start_ind)
#MessageBox.showPlain("%s" %start_time,"Test")
end_time, run_ind = get_datetime_from_sys(end_ind)
#MessageBox.showPlain("%s" %end_time,"Test")

if not os.path.isfile(dss_filename):
	MessageBox.showError("DSS was unable to open the RAS DSS file. \n The filename provided is: %s \n Please check to see that your RAS project is exporting data to this file!"%(dss_filename),"DSS unable to open file!")

dss_file =HecDss.open(dss_filename,start_time,end_time)
try:
	math = dss_file.read(pathname)
	dss_file.done()
	tsc = math.getData()
	values = tsc.values
	times = tsc.times

	time_str = HecTime()
	time_list = []
	for t in times:
		time_str.set(t)
		time_list.append(HecTime.dateAndTime(time_str))

	text_file = open(output_filename, "w")
	text_file.write("Data: \n %s \n" %values)
	text_file.write("Time String: \n %s \n" %time_list)
	text_file.write("Time Minutes: \n %s" %times)
	text_file.close()
except:

	text_file = open(output_filename, "w")
	text_file.write("NO DATA")
	text_file.close()
	MessageBox.showError("DSS was unable to find the path associated with the file. \n The DSS File is: %s \n The Path is: %s \n Please check to see that your DSS file contains the specified path!"%(dss_filename,pathname),"DSS unable to find computed path!")



