# XRB-phot
pipeline for SMARTS xrb optical data

This repository contains code for the prototype X-Ray Binary Photometry Pipeline. Currently the pipeline can reduce optical data, but you might follow along and write a similar program to do the IR data.

An exhaustive discussion on how to use the code that is written, along with examples is avaliable in this repositories wiki. Check it out by clicking the little tab on the right with the book icon on it!

# Dependencies 

For this code to work, you will need to install some other supporting python packages. It is important to ensure the following packages are installed for your /opt/anaconda/bin/python distribution of python. I will show you explicitly how to download these to your system, below.

The required programs are [alipy](http://obswww.unige.ch/~tewes/alipy/index.html)-which itself requires [asciidata](http://www.stecf.org/software/PYTHONtools/astroasciidata/) and [pyfits3](http://www.stsci.edu/institute/software_hardware/pyfits/)-and [pyraf](http://www.stsci.edu/institute/software_hardware/pyraf). 

Other required packages are [astropy](http://www.astropy.org/), [numpy](http://www.numpy.org/), and [pandas](http://pandas.pydata.org/), *although* these ship standard with the anaconda distribution of python, and you ought to have them already.

Before you actually start downloading anything, I would recommend you make sure you do not have these already installed. To check if you have a package installed or not, you can simply start up python and attempt to import it. For example if you can do `import alipy` and you get no error messages, there is no need for you to download alipy, asciidata, or pyfits. 

Alternatively, if you execute the above command and get an error message saying 'ImportError: no module named pyfits' or 'ImportError: no module named asciidata', you likely already have alipy on your system *but not asciidata or pyfits*, respectively. 

Lastly, if you execute the above command and get an error message saying 'ImportError: no module named alipy', you will need to import alipy, asciidata, and pyfits. 

## installing alipy and asciidata

The source code for alipy, and asciidata are in this repository, so I saved you the trouble of googling them and downloading them to your disk. To set up your python with them, follow these instructions.

```shell
cd alipy/
/opt/anaconda/bin/python setup.py install --user #this installs alipy. wait for the output to finish
cd ../asciidata-1.1.1/
/opt/anaconda/bin/python setup.py install --user #this installs asciidata. wait for the output to finish
cd ../
```

## installing pyfits and pyraf

you can use pip to install these packags. Installing packages from pip is **much** easier than installing them from a setup.py file. Do the following in your terminal (you can actually be in whaterver directory you want for this part)

```shell
/opt/anaconda/bin/pip install pyfits --user #this installs pyfits. wait for the output to finish
/opt/anaconda/bin/pip install pyraf --user #this installs pyraf and its dependencies therein. wait for the output to finish.
```

The source code for pyfits and pyraf are also included in this repository for completeness. In principle, it is possible to install them from source following the procedure for alipy and asciidata. However, I *highly* recommend you use pip instead. 
