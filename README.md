AWOT
===============

AWOT is a toolkit of utilities to analyze and visualize weather
observations taken by aircraft.

Authors:

Nick Guy (nick.guy@uwyo.edu)
Andrew Lyons

Created:  11 Sep 2014   Refactored various scripts into a dedicated package

## Installation
The latest source code for AWOT can be obtained from the GitHub repository,
[https://github.com/nguy/AWOT](https://github.com/nguy/AWOT).

Either download and unpack the zip file of the source code or use git to checkout the repository

```python
git clone https://github.com/nguy/AWOT
```
To install in your home directory, use:
```python
python setup.py install --user
```
To install for all users on Unix/Linux:
```python
python setup.py build
sudo python setup.py install
```
 
## Structure
This package is currently divided as follows:

io - Input and output of specific data files
Currently supported:
	NOAA P-3 tail Doppler and lower fuselage radar native radar coordinates.
	NOAA P-3 tail Doppler and lower fuselage radar gridded coordinates, specifically produced by the windysn program.
	NOAA P-3 and UND flight level data.
	LATMOS Falcon NetCDF files which contain both W-band radar and flight information.
	U. Wyoming King Air W-band radar Level 2 NetCDF files.

graph - Produces plots.  Horizontal plots are overlaid on a basemap instance 
and vertical plots are regular 2D plots.

display - Very immature visualization routines.  Work is in progress.

src - Processing software for NOAA P-3 tail radar data.
 
## Usage
Initiate the overall Airborne class

```python
from awot.Airborne import Airborne

fl = Airborne()
```

You can then add data files as such:

```python
fl.get_radar_data(fname='pr_radar_file.nc', platform='p-3', file_format='netcdf', 
                  instrument='tdr_sweep')            
```

```python
fl.get_flight_data(fname='flightfile.nc', platform='p-3', file_format='netcdf')                
```

You should also be able to Initiate the individual data classes:

```python
from awot.io.RadarDataFile import read_radar

radar = read_radar(filename=file, platform='p-3', file_format='netcdf', instrument='tdr_sweep')
```

Data can then be plotted:
```python
from awot.graph.radar_sweep import RadarSweepPlot

rgp = RadarSweepPlot(fl, instrument='tdr_sweep')

rgp.plot_track_relative('dBZ')

from awot.graph.flight_level import FlightLevel

flp = FlightLevel(fl)

flp.plot_trackmap(color_by_altitude=True, track_cmap='spectral', addlegend=True, addtitle=True)
```

## Dependencies

Developed on the Anaconda distribution (1.9.1 & Python 2.7.7).
[Anaconda](https://store.continuum.io/cshop/anaconda/)

It uses a typical scientific python stack:
[Numpy](http://www.scipy.org)
[Scipy](http://www.scipy.org)
[matplotlib](http://matplotlib.org)
[netCDF4](http://code.google.com/p/netcdf4-python)

Development used the 

It also uses the [Py-Art](https://github.com/ARM-DOE/pyart) package for importing some radar data.
And [nappy](https://pypi.python.org/pypi/nappy/1.1.0) to read NASA Ames formatted ASCII files.

## Notes

Future enhancements will include the ability to read data files from other platforms.
  
At present this package is primarily for visualization.  However, work is underway to connect 
many processing routines that allow processing from start to finish with airborne data.

Primary development was with the NOAA P-3 aircraft.  
Ability to display UND Citation and LATMOS Falcon aircraft data was added in, but is not actively developed.
Future development will include U Wyoming King Air flight, radar, and lidar capabilities.

A known issue exists in reading certain King Air radar files.  This is under investigation. 

There is a current dependency on nappy for reading NASA AMES formatted ASCII files. 
This package appears to no longer be supported, though easy to attain.   



