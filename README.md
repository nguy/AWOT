AWOT - Airborne Weather Observations Toolkit
===============

AWOT provides a toolkit of utilities to read and visualize weather observation files collected via aircraft platforms.


##Installation
The latest source code for AWOT can be obtained from the GitHub repository,
[https://github.com/nguy/AWOT](https://github.com/nguy/AWOT).

Either download and unpack the zip file of the source code or use git to checkout the repository

```
git clone https://github.com/nguy/AWOT
```
To install in your home directory, use:

```
python setup.py install --user
```
To install for all users on Unix/Linux:
```
python setup.py build
sudo python setup.py install
```

##Package Structure
_(Documentation to follow soon!)_
###io
Input and output of specific data files<br><br>
In-Situ Flight Data:<br>
&nbsp;&nbsp;* NetCDF<br>
&nbsp;&nbsp;* NASA AMES FFI 1001 format<br>
Remote Sensing:<br>
&nbsp;&nbsp;* Univ. Wyoming Cloud Radar (W-band) Level 2 NetCDF files.<br>
&nbsp;&nbsp;* Univ. of Wyoming Cloud Lidar.<br>
&nbsp;&nbsp;* NOAA P-3/G-IV tail Doppler and lower fuselage radar native radar coordinates.<br>
&nbsp;&nbsp;* NOAA P-3 tail Doppler and lower fuselage radar gridded coordinates produced by the windysn program.<br>
&nbsp;&nbsp;* LATMOS/SAFIRE Falcon NetCDF files which contain both W-band radar and flight information.<br>
&nbsp;&nbsp;* Any surface-based radar file readable with Py-ART.

###graph
Produce horizontal or vertical plots.  Horizontal plots are overlaid on a Basemap instance and vertical plots are regular 2D plots. Experimental 3D plotting is also provided.

###utility
Routines for matching flight track to points in a volume of data, creating radar CFAD diagrams. In addition, a number of helper/conversion routines can be found here.

###src
Processing software for NOAA P-3 tail radar data. Work in progress, not operational yet.

##How-to
Many [examples](https://github.com/nguy/AWOT/tree/master/examples) Jupyter notebooks are provided to familiarize yourself with AWOT.

Read a NOAA P-3 tail radar sweep file:

```
from awot import io
radar = io.read_tdr_sweep(fname='p3_radar_file.nc')
```

Read a flight data file in NetCDF:
```
from awot import io
flight = io.read_netcdf(fname='flightfile.nc', platform='p-3', file_format='netcdf')
```

Read a binary windsyn file:
```
from awot import io

radar = io.read_windsyn_binary(filename=file, platform='p-3', file_format='netcdf', instrument='tdr_sweep')
```

Data can then be plotted:
```python
from awot import graph

rgp = graph.RadarSweepPlot(fl, instrument='tdr_sweep')

rgp.plot_track_relative('dBZ')

from awot import graph

flp = FlightLevel(flight)

flp.plot_trackmap(color_by_altitude=True, track_cmap='spectral', addlegend=True, addtitle=True)
```

##Dependencies

A typical scientific python stack is employed:
    [Numpy](http://www.scipy.org)
    [Scipy](http://www.scipy.org)
    [matplotlib](http://matplotlib.org)
    [netCDF4](http://code.google.com/p/netcdf4-python)

It is highly recommended to use [Anaconda](https://store.continuum.io/cshop/anaconda/) python
distribution that provides easy access to all of these packages along with much more.

Optional:
    [Py-ART](https://github.com/ARM-DOE/pyart) package is used for reading some radar data.
    [Simplekml](http://www.simplekml.com) used to produce KMZ file output.




##Notes
At present this package is primarily for visualization.  However, work is underway to connect
many processing routines that allow processing from start to finish with airborne data.

Development to provide the ability to process Wyoming Cloud radar and lidar files is underway.

AWOT was developed on MacOSX 10.9 - 10.10.5 using Python 2.7, though it should be Python 3 compliant.

##Contributors

Nick Guy (nick.guy@uwyo.edu)

Timothy Lang

Andrew Lyons

