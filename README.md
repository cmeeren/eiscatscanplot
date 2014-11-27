eiscatscanplot
==============

A python script for very easy realtime and offline plotting of EISCAT scans of all types (azimuth, elevation/MSP, and mixed). Currently hard-coded to work on EISCAT Svalbard Radar (ESR), let me know if you need it changed to work for another site (shouldn't take much work).

Usage
-----

Copy `eiscatscanplot.py` to anywhere on your computer, open a terminal, go to the folder and run it using

    `python eiscatscanplot.py`
    
The script will auto-detect the latest 32m folder in the ESR data directory (this can be overridden when you run the script). You may change a few constants at the top of the script to adjust the script behaviour, but this is generally not needed.

The script will update plots in realtime, and when it detects the end of a scan, it will save the figure and start a new plot.

Requirements
------------

* Python 2.7 (tested on Windows 7 64-bit and Ubuntu 14.10 64-bit)
* numpy
* matplotlib
* basemap (from mpl_toolkits)

On a normal Python 2.7 install on e.g. Ubuntu, these can be installed by running

    sudo apt-get install python-numpy python-matplotlib python-mpltoolkits.basemap
