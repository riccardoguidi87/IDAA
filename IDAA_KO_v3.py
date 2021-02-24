#### IDAA KO results

## VERSION 3: this upgrade allows the search of multiple KOs at the same time. It uses as input a table of a defined structure
## where well location of WT and KO are given, as well as expected WT peak size and tollerance range
## the return is a table with results

from __future__ import division
import re
import os
import time
import math 

## FUNCITONS ##
def Find_WT_Peak(list1,wt_peak_size):
#list1 should have the structure of the Dictionary value: |size,height|size,height|size,height|...
#retur the size of the highest peak (assumed to be the WT) and it's relative abundance in the peak range 
#VERSION 2: if the list1 element contains no peaks within expected range, break the conde and give error message
	f = list1.lstrip("|").split("|")
	upper = wt_peak_size + 40
	lower = wt_peak_size - 40
	height = 1
	sum_of_peaks = 1
	highest_peak_size = 1
	for i in f:
		# is the peak in range? if so, store height and store height in SUM of heights
		if  lower < float(i.split(",")[0]) < upper:
			sum_of_peaks += float(i.split(",")[1])
			if height < float(int(float(i.split(",")[1]))):
				height = float(int(float(i.split(",")[1])))
				highest_peak_size = float(i.split(",")[0])
	#if float(height) == 1:
		#print "I didn't find a peak of the expected size in the WT reference well" 
	relative_abundance = float(height) / float(sum_of_peaks)
	return highest_peak_size , round(relative_abundance, 3)

def KO_efficiency(list1,wt_peak_size):
#list1 should have the structure of the Dictionary value: |size,height|size,height|size,height|...
#wt_peak_size is the value obtained by Find_WT_Peak funciton
#return the relative abundance of the WT peak in the area surrounding the peak. This funciton requires 
#the identificaiton of the exact WT peak based on the WT run - see Find_WT_peak function
# note that the peak-range here is more restrictive then in the WT search, as we expect CRISPR cuts to be below 10bp

	f = list1.strip("|").split("|")
	upper = wt_peak_size + 10
	lower = wt_peak_size - 10
	sum_of_peaks = 1
	size = 1
	delta_in_hold = 100
	height = 1
	possible_wt_size = 1
	for i in f:
		# is the peak in range? if so, store height and store height in SUM of heights
		if  lower < float(i.split(",")[0]) < upper:
			sum_of_peaks += float(i.split(",")[1])
			delta = float(i.split(",")[0]) - float(wt_peak_size)
			delta = abs(delta)
			if delta < delta_in_hold: # here we test if the peak is closer to the WT size than a previously peak
				delta_in_hold = delta
				possible_wt_size = float(i.split(",")[0])
				height = float(i.split(",")[1])
	#if float(height) == 1:
		#print "although I found the WT peak in the WT sample, I didn't find the same peak in the KO sample, so I cannot calculate the relative KO efficienty"
	relative_abundance = float(height) / float(sum_of_peaks)
	return possible_wt_size, round(relative_abundance, 3)

# get info on how data are organized in the original Thermo table
with open('Table1.csv') as f: 
	line = f.readline().rstrip("\r\n").split(",")
	position_SampleName = line.index('"Sample File Name"')
	position_DyeColour = line .index('"Dye Color"')
	position_Size = line.index('"Size"')
	Position_Height = line.index('"Height"')

## create empty dictionary
Peaks_Dict = {}

## pre-populate with all the Sample name files. This way you will create unique Key entries, all empty
## VERSION 3: I pre-populate the Peaks_Dict dictionary with a fake "empty" peak of value 1,1 to bypass an error that
## happens when using Functions if the dictionary value is empty
Table1 = open("Table1.csv")
next(Table1) # this comand skip the first line of the file
for line in Table1:
	SampleName = line.rstrip("\r\n").replace('"',"").split(",")[0]
	Peaks_Dict[SampleName] = "|1,1"
Table1.close()

## populate Dicrionary Peaks_Dict directly from the orignal Thermo table
Table1 = open ("Table1.csv")
next(Table1) # this comand skip the first line of the file
for line in Table1:
	linelist = line.rstrip("\r\n").replace('"',"").split(",")
	#print repr(line)
	if linelist[position_DyeColour] == 'BLUE':
		if float((linelist[position_Size])) >= 100:
			if float((linelist[Position_Height])) >= 100:
				SampleName = linelist[0]
				Peaks_Dict[SampleName] = Peaks_Dict[SampleName] + "|" + str(linelist[position_Size]) + "," + str(linelist[Position_Height].rstrip('\n'))
Table1.close()

## at this stage, you should have Dictionary that contains one key per well. Each Key has a value that correspond
## to the peak size and height (anything larger then 100bp and highr then 100) 
## Key = | size, height| size, height| size, height| size, height 

#### PART 2: Interrogate the Dictionary via iterations
## use a User_input table from which outsourse all the necessary info for the search

output = open("output.csv","w")
output.write("PlateID,WT_LOCATION,KO_LOCATION,PEAK_SIZE,Tollerance Range,Exp,crRNA,WT_peak,WT_peak_abd,WT_peak_inKOcells,WT_peak_abd_inKOcells\n")

User_input = open("User_Input.csv")
next(User_input)
for line in User_input:
	PlateID = line.split(",")[0]
	WT_well = line.split(",")[1]
	KO_well = line.split(",")[2]
	WT_PEAK_SIZE = float(line.split(",")[3])
	Tollerance = float(line.split(",")[4])

	#run loop
	for key, value in Peaks_Dict.items():
		if key.split("_")[2] == PlateID:         ## CAREFUL HERE check SampleNames format!!
			if key.split("_")[4] == WT_well:         ## CAREFUL HERE check SampleNames format!!
 				WT_peak = Find_WT_Peak(Peaks_Dict[key],WT_PEAK_SIZE)[0]
				WT_peak_abd = Find_WT_Peak(Peaks_Dict[key],WT_PEAK_SIZE)[1]
	
	for key, value in Peaks_Dict.items():
		if key.split("_")[4] == KO_well:   # these conditions was removed: "WT_peak != 1 and" and also: "key.split("_")[2] == PlateID"
			WT_peak_inKOcells = KO_efficiency(Peaks_Dict[key],WT_peak)[0]
			WT_peak_abd_inKOcells = KO_efficiency(Peaks_Dict[key],WT_peak)[1]
			output.write(line.rstrip('\r\n') + "," + str(WT_peak) + "," + str(WT_peak_abd) + "," + str(WT_peak_inKOcells) + "," + str(WT_peak_abd_inKOcells) + '\n')

User_input.close()
output.close()


