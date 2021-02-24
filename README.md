## Automated Analysis of CRISPR KO efficiency via Indel Detection by Amplicon Analysis (IDAA)
### A Python script that compares peaks of IDAA electrophoresis analysis to automatically infer KO frequency of CRISRP treatement   

This program takes two tables as input:
1. **ThermoFisher Cloud Electrophoresis Tabulation** (TFET): you can obtain such table when you run your electrophoresis files (.fsa) using a [Cloud Thermo](https://www.thermofisher.com/order/catalog/product/A26811?SID=srch-srp-A26811#/A26811?SID=srch-srp-A26811) account, or any registered Thermo Electrophoresis softwares that can process .fsa file formats (that is the raw data from electrophoresis)
2. **User-Input-Table**: containing the necessary instructions for the identification of a specific peak, or set of peaks, and its quantification

Both tables are *.csv*

In the basic version of the program, the code searches for one peak, indicated in the User-Input-Table, inside the TFET. 

This program assumes you are running a 96WP-format of electrophoresis. You don't need to have a TFET containing all 96 sample, 
but each samples must be labelled with a row (A-H) and column system (1-12) that corresponds to a 96WP.

In the .fsa file, metadata contain the name of each sample, which will be passed to the TFET to label each individual sample. In this  version of the program, the nomenclature of each sample **must follow a precise order** (see below). The nomenclature I use for all my electrophoresis runs - and therefore the one this program can handle - looks like this:

```SampleName``` : "RUN200225B_A01_A1_015_extrainfo.fsa"

Where we have in order:

```SampleName``` :  *plateID _ wellcode _wellcode_samplenumber_extrainfo.fsa*

I am not sure why, but TFET contains a repetition of the ```wellcode``` (row and Col values combined). The program needs these two ```wellcode``` to be present in the ```SampleName``` in this precise order.

If the ```SampleName``` is prepared differently, and cannot be easily changed, take a look at ```#run loop``` part of the Python script to fix the problem: 

```python
# run loop
	for key, value in Peaks_Dict.items():
		if key.split("_")[0] == PlateID and key.split("_")[2] in wells_range:  ## CAREFUL HERE check SampleNames format!!
		#if key.split("_")[3] in wells_range:  ## USE THIS IF you have only one plate
			WellID = key.split("_")[2]
```

... the location ```key.split("_")[2]``` of the split string refers to the location in ```SampleName``` where to find the ```wellcode``` in format ND (see above). If ```wellcode``` is present in another location, adjust accordingly. Also, if the split-value isn't "_", that can be changed here. Also ```key.split("__")[4]``` make sure this location in the SampleName hold the well ID in a format that is ‘A1’, not ‘A01’

WellID = key.split("_")[4]
make sure this IS the location of the well ID in SampleName a format that is ‘A1’, not ‘A01’

- The program does take into consideration you may have more then one plate in Table1.csv, and each plate may contain more then one PP. Make sure that in User_Input.csv you have correctly placed the plate ID (it does’t have to be a number, it can also be a string - most recent agreement in plate labeling from Dixon, the PlatesID are continuous letters A,B,C…. etc… ).

### Instructions to Run Analysis

- A copy of the original TFET and the User-Input-Table must be present inside the same folder where the Python program runs
- Rename TFET as: “Table1.csv” - be exact with the name!
- The User-Input-Table must be named as "User_Input.csv", and should have this outlook:

### User-Input-Table Outlook


| PLATE      | WT_LOCATION  | KO_LOCATION| PEAK_SIZE | Tollerance Range | Exp   | crRNA     |
| ---------- | ------------ | ---------- | --------- | ---------------- | ----- | --------- |
| RUN200916A | A1           | A2         | 300       | 20               | RG222 | gRNA.21   |
| RUN200916B | A1 	    | A3         | 300       | 20               | RG222 | gRNA.22   |
| RUN200916C | A1   	    | A4         | 300       | 20               | RG222 | gRNA.23   |

Notes on the User Input Table:

- Even if you have run only one electrophoresis plate (in the example above I run three), the name of plate must always be present in User Input Table, as well as in the ```SampleName``` of the TFET

- When modify the User_Input_Table using Excel, make sure you ERASE ANYTHING (even a blank cell) from underneath the last row (the program sometimes see the leftovers and tries to continue reading rows). Save table as *.csv*
- The program takes into consideration continues and ordered use of wells in the 96WP sent for electrophoresis - the program may still be working fine if some wells are to be missing from the run (as in, Table1.csv does not contains a raw for each of the 96 wells). 
- Beware of the way the facility lab that runs the IDAA plate labels the samples (usually first column in the exported Thermo Table): In the #run loop section, the program relies on a specific order of information present in the “Sample File Name” column of Table1.csv.


### Explanations of output table

| PlateID | WT\_LOCATION | KO\_LOCATION | PEAK\_SIZE | Tollerance Range | Exp     | gRNA | WT\_peak | WT\_peak\_abd | WT\_peak\_inKOcells | WT\_peak\_abd\_inKOcells |
| plate2 | B2 | B3 | 319 | 20 | RG387.5 | gR.23 | 316.91 | 0.466 | 316.95 | 0.298 |
| ------ | -- | -- | --- | -- | ------- | ----- | ------ | ----- | ------ | ----- |
| plate2 | B2 | C3 | 319 | 20 | RG387.5 | gR.23 | 316.91 | 0.466 | 316.95 | 0.305 |
| plate2 | B2 | D3 | 319 | 20 | RG387.5 | gR.23 | 316.91 | 0.466 | 317    | 0.31  |
| plate2 | B2 | E3 | 319 | 20 | RG387.5 | gR.23 | 316.91 | 0.466 | 316.95 | 0.308 |
| plate2 | B2 | F3 | 319 | 20 | RG387.5 | gR.23 | 316.91 | 0.466 | 316.95 | 0.3   |
| plate2 | B2 | G3 | 319 | 20 | RG387.5 | gR.23 | 316.91 | 0.466 | 316.86 | 0.29  |

#Same Values foound in the **User-Input-Table**:
PlateID	= name of plate run. The software can analyze multiple plates at once, as long as the plate number is includede within ```SampleName``` of TFET
WT_LOCATION = this is the well coordinate where to expect a peak corresponding to a WT allele. This is usually the location of primer pairs that amplyfy the exon of interest from a CRISPR-CTR sample.
KO_LOCATION = this is the well coordinate where a CRISPR should have cut a gene. The primer pair used in this well should be identical to the primer-pair used int to amplify the peak in "WT_LOCATION".
PEAK_SIZE = What size do you expect your peak to be at?
Tollerance Range = the error range within which we can expect to see the WT peak
Exp = personal experiment number (not relevant to the program)
gRNA = name of the gRNA used (not relevant to the program)

#Output Values:
WT_peak	= this is the size of the WT allele peak identify in the CTR sample (see WT_LOCATION)
WT_peak_inKOcells = this is the size of the WT allele peak identify in the CRISPR sample (see KO_LOCATION). **WT_peak_inKOcells** and **WT_peak** should be identical, indicating that both peaks were correctly identified by the program. 
WT_peak_abd_inKOcells = This value is calculated exclusivelly using peaks detected in the CRISPR-KO. In the CRISPR-KO well, the WT peak is first identified based on the size found in **WT_peak**. Then, additional, surrounding peaks are identified in the CRISPR-KO well. This value is the ration between the height of the WT peak, vs all other non-WT peaks within 10bp of range.

















