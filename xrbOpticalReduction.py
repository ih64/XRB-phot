import pandas as pd 
import numpy as np
import alipy, glob, os, stat, subprocess
from astropy.io import ascii, fits
from pyraf import iraf
#do this stuff down here so IRAF doesn't save any parameters, and to make important tasks avaliable 
iraf.prcacheOff()
iraf.noao.obsutil()
iraf.imred()
iraf.ccdred()

def alignProcessed(target,fltid,ref_image,outdir):
	'''align the processed ccd images using alipy, a python module I found that calls
		source extractor and iraf imalign tasks 
	INPUT:
	target: a string indicating the target you want to reduced
	fltid: the andicam filter for which these data were taken in.
	ref_image: a string indicating the path to the refence image used to align the other images
	outdir: a string indicating the path to save the new aligned fits images to
	OUTPUT:
	new fits images that are aligned to the reference images. 'geregister' is depended 
	to the file name of each input image'''
	#grab the paths to the images we are going to align
	images_to_align =sorted(glob.glob("/net/xrb-archive/usb-data/"+target+"/fitsimages/"+fltid+"/*.fits"))

	#everything below here i copied from the alipy demo http://obswww.unige.ch/~tewes/alipy/tutorial.html
	#this line uses source extractor to identify stars in the images
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	# That's it !
	# Put visu=True to get visualizations in form of png files (nice but much slower)
	# On multi-extension data, you will want to specify the hdu (see API doc).

	# The output is a list of Identification objects, which contain the transforms :
	# These tell you how rotate and translate the images so they are aligned to the referene image
	for id in identifications: # list of the same length as images_to_align.
		if id.ok == True: # i.e., if it worked
			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
			# id.trans is a alipy.star.SimpleTransform object. Instead of printing it out as a string,
			# you can directly access its parameters :
			#print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
			#print id.trans.matrixform()
			#print id.trans.inverse() # this returns a new SimpleTransform object
		else:
			print "%20s : no transformation found !" % (id.ukn.name)

	# Minimal example of how to align images :
	outputshape = alipy.align.shape(ref_image)
	# This is simply a tuple (width, height)... you could specify any other shape.

	#finally, for each image where a transform was found, create a new image where the data are transformed
	#to be alinged with the reference image, and print it to 'outdir', the last argument of this function
	for id in identifications:
		if id.ok == True:
			# Variant 2, using geomap/gregister, correcting also for distortions :
			alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars,verbose=False,
				shape=outputshape, makepng=False, outdir=outdir)
	return

def avgFWHM(image,coords):
	'''use the iraf task psfmeasure to compute the average fwhm of a few stars
	INPUT: 
	image: the fits image you want to measure the fwhm for
	coords: a text file listing the x,y pixel coordinates of stars in image to use to compute the fwhm.
			I'd recommend choosing ~10, bright stars, evenly sampled across the frame.
	OUTPUT: a string of the average fwhm
		also your graphics terminal will open up with a gross iraf graphics terminal
		If you are interested you can look at it to see info like elipticity and average fwhm of the image.
		
	If you run this in batch mode and consistently are seeing fwhm values >~ 5, it might mean you need to choose different
	stars in your coords file, or check that the coordinates are accurate.
	'''

	#we need to create a text file that has the letter 'q' in it
	#this gets fed to psfmeasure as the graphcur input. it forces iraf to quit out of the terminal
	#this way you can compute the fwhm with out the user having to click on anything
	with open('graphcur','w') as f:
		f.write('q')

	#use the iraf task psfmeasure to get the fwhm
	#the output is not machine readable friendly, but basically we are after the last thing it prints out
	#this is the avg fwhm measured
	stdOutput=iraf.psfmeasure(image,imagecur=coords,display='no', graphcur='graphcur',Stdout=1, logfile='')
	FWHM=stdOutput[-1].split()[-1]

	#clean up the graphcur textfile we made
	iraf.delete('graphcur')

	return FWHM

def prepDAOPHOT(fwhm):
	'''create the input tables daophot needs to do its job
	you need a priori info about the fwhm and detector characteristics to use daophot well
	this program figures that out and prepares the tables for you to do psf photometry on
	INPUT:
	fwhm: the full width half max of an image
	OUTPUT:
	daophot.opt: a text file that lists the parameters for daophot
	photo.opt: a text file that lits the parameters daophot uses to do aperture photometry
	allstar.opt: a text file that allstar uses to do psf photometry'''

	with open('daophot.opt','w') as dao:
		#use the fwhm to determine the fwhm, fit, and psf
		dao.write('FWHM='+str(fwhm)+'\n')
		dao.write('FIT='+str(fwhm)+'\n')
		dao.write('PSF='+str(float(fwhm)*3.0)+'\n')
		#put in andicam characteristics
		dao.write('READ=6.5\n')
		dao.write('GAIN=2.3\n')
		#set the threshold to 3.5 sigma
		dao.write('TH=3.5\n')
		#use gausian analytic model
		dao.write('AN=1\n')
		#set the low signal as 10 sigma so we dont read into the noise
		dao.write('LOWBAD=10\n')
		#andicam is linear up to this many adu
		dao.write('HIBAD=45000\n')
		#don't ask for user input
		dao.write('WATCH=-2.0\n')
		#the psf should remain the same across the detector
		dao.write('VAR=0\n')

	with open('photo.opt','w') as photo:
		#set the photometry radius to the fwhm
		photo.write('A1='+str(fwhm)+'\n')
		#inner sky radius to 3*fwhm
		photo.write('IS='+str(float(fwhm)*3.0)+'\n')
		#outer sky radius to 4*fwhm
		photo.write('OS='+str(float(fwhm)*4.0)+'\n')

	with open('allstar.opt','w') as allstar:
		#use the fwhm as the psf fitting radius
		allstar.write('fit='+str(fwhm)+'\n')
		#inner skyr adius to 3*fwhm
		allstar.write('isky='+str(float(fwhm)*3.0)+'\n')
		#outer sky radius to 4+fwhm
		allstar.write('osky='+str(float(fwhm)*4.0)+'\n')
		#dont ask for user input
		allstar.write('watch=0\n')
		#allow all star to redetermine centroids to give a better quality of fit
		allstar.write('redet=1\n')

	return

def daophotWrap(filename):
	'''creates a shell script that will call daophot to run find, pick, and psf 
	on the image "filename"
	INPUT:
	filename: a string of the filename you want to run daophot on
	OUTPUT:
	daophotGo.sh: a shell script that wraps around daophot
	'''
	with open('daophotGo.sh','w') as f:
		f.write('daophot <<__DAOPHOT-END__\n')
		f.write('attatch '+filename+'\n')
		f.write('find\n')
		f.write('1,1\n')
		f.write(filename+'.coo\n')
		f.write('y\n')
		f.write('phot\n')
		f.write('photo.opt\n')
		f.write('\n')
		f.write(filename+'.coo\n')
		f.write(filename+'.ap\n')
		f.write('pick\n')
		f.write(filename+'.ap\n')
		#pick 20 stars, magnitude limit 20
		f.write('20,20\n')
		f.write(filename+'.lst\n')
		f.write('psf\n')
		f.write(filename+'.ap\n')
		f.write(filename+'.lst\n')
		f.write(filename+'.psf\n')
		f.write('exit\n')
		f.write('__DAOPHOT-END__')
		f.write('\n')

	#give executable permission to the shell script
	os.chmod('daophotGo.sh',0755)
	return

def allstarWrap(filename):
	'''creates a shells cript that will call allstar to do psf photometry on the image "filename"
	INPUT:
	filename: a string of the filename you want to do psfphotometry on
	OUTPUT:
	allstarGo.sh: a shell script that wraps around allstar
	'''
	with open('allstarGo.sh','w') as f:
		f.write('allstar <<__ALLSTAR-END__\n')
		f.write('\n')
		f.write(filename+'\n')
		f.write(filename+'.psf\n')
		f.write(filename+'.ap\n')
		f.write(filename+'.als\n')
		f.write('\n')
		f.write('__ALLSTAR-END__\n')
		f.write('\n')

	#give executable permissions to shell script
	os.chmod('allstarGo.sh', 0777)
	pass

def findSources(alsfile,coords,tol):
	'''look in the output of allstars, alsfile, and find the entries of particular stars,
	whose x,y pixel coordinate positions are given in the textfile, coords. 
	Uses a vectorized nearest neighbor search
	this function was made to be called by others like makePSFPhotDF below, to help return a photometry table for a target
	by itself, this function only returns photometry information for one target
	INPUT:
	alsfile: the file name of the allstar output. typically this has a .als multi-extension
	coords: the file name of a a text file that lists the coordinates of the stars you are interested in
	tol: tolerance distance in pixels. If we cannot find any stars with distances below the tol
	we don't remember count its daophot data and instead return nans
	OUTPUT:
	photBios: a python dictionary with photometry info for these sources 
	'''
	#first read in the psf photometry output file to a pandas dataframe
	alsPhot=pd.read_table(alsfile, skiprows=[0,1,2],
		names=['ID','x','y','mag','err','sky','nit','chi','sharp'], sep='\s+')

	#now read the coords file into a dataframe
	coords=pd.read_table(coords, names=['x','y','junk'], sep='\s+')

	#count up how many sources are in the coords list
	numSources=len(coords)

	#prepare a numpy array to hold the dataframe indicies for the stars in the coords list
	starIlocs=np.zeros(numSources)

	#loop over the stars in the coords dataframe. for every entry, find its nearest neighbor
	#in the alsPhot dataframe
	for i in xrange(numSources):
		coordsView=coords.iloc[i]
		x=coordsView.x
		y=coordsView.y
		#calculate the distances between this star in coordsView and all the stars in the als dataframe
		distances=np.sqrt(np.square(alsPhot.x - x) + np.square(alsPhot.y - y))
		#if the smallest distance is within the tolerance range, keep this star's index number in 
		#the als data frame in the starIloc array
		if distances.min() < tol:
			starIlocs[i]=np.sqrt(np.square(alsPhot.x - x) + np.square(alsPhot.y - y)).argmin()
		#otherwise just put a nan there
		else:
			starIlocs[i]=np.nan

	#prepare a python dict to keep the data and meta data from the als file on the stars we found
	photBios={}

	#loop over the stars we found, look them up in the alsPhot dataframe, and save their info 
	#in the photBios dict.
	for i in xrange(numSources):
		#if we were able to find a star with distance within tol save its alsPhot data in photBios
		if np.isfinite(starIlocs[i]):
			alsPhotView=alsPhot.iloc[starIlocs[i]]
			for j in alsPhotView.index:
				photBios['s'+str(i)+'_'+j]=alsPhotView[j]
		#otherwise just fill everything with nans
		else:
			for j in alsPhot.columns:
				photBios['s'+str(i)+'_'+j]=np.nan

	return photBios

def makePSFPhotDF(filelist,fwhmthresh=8.0,psfcoords='psfcoords.lis',photcoords='photcoords.lis',pickle=True,csv=True):
	'''
	run daopht, allstar, and extract photometry for specific targets, and return a photometry log for specefied targets
	to use it, first cd to /net/xrb-archive/usb-data/TARGET/fitsimages/flt_align/ where flt is the filter (either B,V,R, or I)
	INPUT:
	filelist:a python string of the aligned images you want to do psf photometry of
	fwhmthresh: a float for the pixel value of the largest allowed fwhm. if a frame in filelist has a fwhm > than fwhmthresh,
		daophot and allstar will not be run on the image. The default value of 8.0 pixels corresponds to about 3 arc seconds of seeing.
	psfcoords: a text file with two columns, specifying x and y pixel coordiates for a few stars to use to measure the fwhm
		in images in filelist the program will iterate over
	photcoords: a text file with two columns, specifying x and y pixel coordinates for stars whose psf photometry will be extracted
		from the *.asl file that is created by allstar for the images in file list that the program will iterate over
	OUTPUT:
	photdf: pandas dataframe that contains header info for the files in filelist like obsdate filename etc,
		as well as all the the info from the *.als for stars specefied by phorcoords
	photdf.pkl: if pickle=True, photdf is pickled and saved to disk in the current directory as 'photdf.pkl'
	photdf.csv: if csv=True, photdf is saved to csv file to disk in the current directory as 'photdf.csv'

	EXAMPLE:
	cd /net/xrb-archive/usb-data/GX339-4/fitsimages/I_align
	filelist=sorted(glob.glob('*.fits'))
	makePSFPhotDF(filelist)
	'''
	row_list=[]

	for f in filelist:
		#daophot won't work on any filenames that have more than 30 characters.
		#if the file name is too long, just skip it completly
		if len(f) < 31:
			hdulist=fits.open(f)
			header=hdulist[0].header
			hdulist.close()
			#get some useful info about the image from its header
			#occasionally the header info is corrupted, so we add error handling
			#just put in np.nan if there is a key error
			try:
				JD=header['JD']
			except KeyError:
				JD=np.nan
			try:
				airmass=header['SECZ']
			except KeyError:
				airmass=np.nan
			try:
				filename=header['FILENAME']
			except KeyError:
				filename=np.nan
			try:
				obsdate=header['DATE-OBS']
			except KeyError:
				obsdate=np.nan
			try:
				exptime=header['EXPTIME']
			except KeyError:
				exptime=np.nan

			#use the psfcoords.lis file to measure the fwhm of f
			fwhm=avgFWHM(f,psfcoords)

			#if the fwhm is lower than the fwhmthresh, proceed with the reduction
			if float(fwhm) < fwhmthresh and float(fwhm) > 1.0:
				#create the .opt files daophot and allstars use
				prepDAOPHOT(fwhm)
				#run daophot on this image
				daophotWrap(f)
				os.system('./daophotGo.sh')

				#sometimes daophot cannot create a psf.
				#in this case, it creates an empty .psf file
				#we must check the length of the .psf file here, and proceed if it succeeded
				psf_line_num=int(subprocess.check_output(['wc', '-l', f+'.psf']).split()[0])

				if psf_line_num > 0:
					#run allstar on this image
					allstarWrap(f)
					os.system('./allstarGo.sh')
					
					#all star might have encoutered a problem. make sure the .als file exists
					#before we try to read into it
					if len(glob.glob(f+'.als')) == 1:
						#get the psfphotoemtry data of the stars in the photcoords.lis file
						photBios=findSources(f+'.als',photcoords,float(fwhm))
						#add the header info we found earlier to the photBios dict
						#put in error handling, sometimes the fits headers will have bad values that cant be converted to floats
						try:
							photBios['JD']=float(JD)
						except ValueError:
							photBios['JD']=np.nan
						try:
							photBios['airmass']=float(airmass)
						except ValueError:
							photBios['airmass']=np.nan
						try:
							photBios['fwhm']=float(fwhm)
						except ValueError:
							photBios['fwhm']=np.nan
						try:
							photBios['exptime']=float(exptime)
						except ValueError:
							photBios['exptime']=np.nan
						photBios['filename']=filename
						photBios['obsdate']=obsdate
						photBios['align_filename']=f
						#add this dict to the running list
						row_list.append(photBios)
					else:
						print "could not find an .als file. "+f+" wont be added to photometry log"

				#if the psf file is empty print a message to the stdout saying we're skippnig psfphotometry
				else:
					print "there is no psf for "+str(f)
					print "skipping psf photometry for "+str(f)

			#if the fwhm was too big, skip the photometry steps and print out an error message
			else:
				print "the fwhm for "+str(f)+" is "+str(fwhm)+" which is bigger than the threshold "+str(fwhmthresh)
				print "skipping photometry for "+str(f)
				
		else:
			print "the file "+f+" has "+str(len(f))+" characters"
			print "daophot will only accept files with 30 characters or less"
			print f+" will not be reduced"

		#clean up the intermediate files before going on to the next file
		flag=os.system('rm *.opt *.sh *.ap *.coo *.lst *.nei *s.fits')
	#shove everything ito a pandas data frame and return it	
	photdf=pd.DataFrame(row_list)
	photdf.to_csv('photdf.csv')
	return photdf