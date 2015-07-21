'''
this is a script to create zeroPoints.csv, a comma seperated value file that contains the different zeroPoints suzane calculated
for dates between 2006 and 2012.

to execute this script, first create a session of ipython.
then change directories to /net/xrb-archive/usb-data/REDUCTION/
then type "run xrbOpticalCalibration"
it will read the zeroPoints tables suzane created (whcih have been copied to /net/xrb-archive/usb-data/REDUCTION/zeroPoints/)
and output /net/xrb-archive/usb-data/REDUCTION/zeroPoints/zeroPoints.csv

you can use this to calibrate the true magnitudes of calibration stars
'''
import pandas as pd
import numpy as np
import glob

suzpts=sorted(glob.glob('/net/xrb-archive/usb-data/REDUCTION/zeroPoints/RESULTS.200[6-9]_*')+glob.glob('/net/xrb-archive/usb-data/REDUCTION/zeroPoints/RESULTS.201[0-1]_*') + glob.glob('/net/xrb-archive/usb-data/REDUCTION/zeroPoints/RESULTS.2012_0[0-6]'))
rowdict={'date':[], 'star':[], 'b':[], 'v':[],'r':[],'i':[], 'cb':[],'cv':[],'cr':[],'ci':[],'xb':[],'xv':[],'xr':[],'xi':[]}
headerdelim='======================================================================================================================='
for month in suzpts:
	#they're not super machine-readable-friendly, we have to take them apart and yank out what we need,
	with open(month,'r') as f:
		text=f.read()
	#header info is in the first two sections above the delimiter,
	data=text.split(headerdelim)[2]
	#each date is broken up by a long string of hyphens,
	dates=data.split('-----------------------------------------------------------------------------------------------------------------------')
	#go up to the second to last element in the list dates, the last one is just an empty string,
	for day in dates[0:-2]:
		words=day.strip().split()
		if words[0][0:2]=='20':
			rowdict['date'].append(float(words[0][:]))
			rowdict['star'].append(words[1])
			rowdict['b'].append(words[2])
			rowdict['v'].append(words[3])
			rowdict['r'].append(words[4])
			rowdict['i'].append(words[5])
			rowdict['xb'].append(words[6])
			rowdict['xv'].append(words[7])
			rowdict['xr'].append(words[8])
			rowdict['xi'].append(words[9])
			rowdict['cb'].append(words[10])
			rowdict['cv'].append(words[11])
			rowdict['cr'].append(words[12])
			rowdict['ci'].append(words[13])

#shove all the data into a pandas dataframe
table=pd.DataFrame(rowdict)
#force the datatype for the dates to be ints, they were made from strings to floats above
table['date']=table['date'].values.flatten().astype(int)
#the mirror was cleaned on 110921, this note sneaks in and screws up the data frame
table=table.drop(table[table['date']==110921].index)
#suzane filled in the string 'na' if the data was not determined. change these to np NaNs in our dataframe
table=table.replace('na',np.nan)
#sometimes the string '---' was used if the data were not determined, swap these out for nans too
table=table.replace('---',np.nan)
table=table.replace('----',np.nan)
#finally, sometimes 0.0 were used if the data were not determined
table=table.replace(0.0,np.nan)
#foce the datatype for the values to floats, they are currently strings.
table[['b','v','r','i','xb','xv','xr','xi','cb','cv','cr','ci']]=table[['b','v','r','i','xb','xv','xr','xi','cb','cv','cr','ci']].astype(float)
table.to_csv('/net/xrb-archive/usb-data/REDUCTION/zeroPoints/zeroPoints.csv', index=False)