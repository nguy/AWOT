"""
awot.graph.plot_tdr_swp
=========================

A group of scripts create various plots of data collected by the NOAA P-3 tail Doppler radar. 

Created by Nick Guy.

"""
# HISTORY::
#  6 Mar 2014 - Nick Guy NOAA/NSSL/WRDD, NRC
# FUNCTIONS::
# polar_sweep - Plot polar coordinate data on polar coordinate axis
# polar_sweep_grid - Plotting transformed data to Cartesian output
# sweep_to_Cart - Polar coordinates transformed to Cartesian
# sweep_aircraft_relative - Polar coord data transformed to aircraft-relative frame
# sweep_track_relative - Polar coord data transformed to track-relative frame
# sweep_earth_relative - Polar coord data transformed to earth-relative frame
#-------------------------------------------------------------------
# Load the needed packages
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import ticker
import numpy as np

import general.library as gl
import general.gplot as gp

# Define various constants that may be used for calculations
#===============================================================
# BEGIN FUNCTIONS
#**===============================================================
def polar_sweep(Var,rot,range,nlevs=30,
               vmin=None,vmax=None,cmap=None,mask_outside=True,
               rng_rings=None,rot_angs=None,
               title=None,cb_flag=True,cb_orient='horizontal',cb_lab=None):
    """Plot a sweep of native (polar) coordinate radar data on polar format plot
 INPUT::
  Var             = Data values to plot
  range           = Range along ray
  rot             = Rotation angle with respect to instrument [degrees]

       THESE NEXT VALUES CUSTOMIZE THE PLOT
  nlevs           = Number of contour levels to plot
  vmin            = Minimum value to display
  vmax            = Maximum value to display
  cmap            = Matplotlib colormap name
  mask_outside    = Boolean flag to mask values outside vmin/vmax
  title           = Title to label plot with, if none given it is omitted
  cb_flag         = True turns on colorbar, False no colorbar
  cb_lab          = Colorbar label, None will use a default label generated from the
                     field information.
  cb_orient       = Orientation of colorbar (vertical or horizontal)
 OUTPUT::
  fig             = Figure to create plot to
  ax              = Polar axis created
  p               = Output plot
 USAGE::
  p, fig, ax = DPJgrid_Horiz(fig,Z,W,Lon,Lat,Ht,TLon,TLat,TAlt,[Zcoord],[Ucoord],[WindVec],[Track],
                     [cminmax],[clevs],[vmin],[vmax],[dlat],[dlon],[proj],
                     [cmap],[title],[pName],[pType],[figsize])
 NOTES::
 Defaults are established during DYNAMO project analysis
    """
# HISTORY::
#  15 Apr 2014 - Nick Guy NOAA/NSSL/WRDD, NRC 
#--------------------------------------------------------
    # Plot the polar coordinate radar data    
    fig, ax, p = gp.plot_polar_contour(Var,rot,range,nlevs=nlevs,
                                      vmin=vmin,vmax=vmax,cmap=cmap,mask_outside=True)

    # Set the range and turn grid on
    ax.set_rmax(1.05 * range.max())
    ax.grid(True)
               
    # Set the title
    if title == None:
        pass
    else:
        ax.set_title(title)
        
    # Set the range rings if set
    if rng_rings == None:
        pass
    else:
       plt.rgrids(rng_rings)
        
    # Set the rotation angles (theta) if set
    if rot_angs == None:
        pass
    else:
       plt.thetagrids(rot_angs)
       
    # Set the colorbar if desired
    if cb_flag:
        cb = plt.colorbar(mappable=p,orientation=cb_orient)
        if cb_lab == None:
            pass
        else:
            cb.set_label(cb_lab)
    
    return fig, ax, p
#**===============================================================
def plot_sweep_grid(Xcoord,Ycoord,Values,ax=None,title=None,
               vmin=-24.,vmax=64.,cmap='jet',grid_on=True,
               xlims=None,ylims=None,xlab=None,ylab=None,
               cb_flag=True,cb_orient='horizontal',cb_lab=None):
    """Plot a sweep of native (polar) projected on a Cartesian plane.
 INPUT::
  X               = Variable to plot along the X-axis
  Y               = Variable to plot along the Y-axis
  Values          = Data values to plot (Ydims,Xdims)

       THESE NEXT VALUES CUSTOMIZE THE PLOT
  ax              = Axis instance on which to be plotted
  title           = Title to label plot with, if none given it is omitted
  vmin            = Minimum value to display
  vmax            = Maximum value to display
  cmap            = Matplotlib colormap name
  xlims           = X-axis limits (min,max)
  ylims           = Y-axis limits (min,max)
  xlab            = X-axis label
  ylab            = Y-axis label
  cb_flag         = True turns on colorbar, False no colorbar
  cb_lab          = Colorbar label, None will use a default label generated from the
                     field information.
  cb_orient       = Orientation of colorbar (vertical or horizontal)
 OUTPUT::
  fig             = Figure to create plot to
  p               = Output plot
 USAGE::
  p = plot_sweep_grid(X,Y,Val,[**args])
 NOTES::
 Plotting convention is to project the data onto a 2D surface looking from the back of
   aircraft forward
    """
# HISTORY::
#  15 Apr 2014 - Nick Guy NOAA/NSSL/WRDD, NRC 
#--------------------------------------------------------
    # Plot the data
    p = ax.pcolormesh(Xcoord,Ycoord,Values,cmap=cmap,vmin=vmin,vmax=vmax)
               
    # Set the title
    if title == None:
        pass
    else:
        ax.set_title(title)
    # Set the axes limits
    if xlims == None:
        pass
    else:
        ax.set_xlim(xlims)
    if ylims == None:
        pass
    else:
        ax.set_ylim(ylims)

        
    # Set the axes labels
    if xlab == None:
        pass
    else:
        ax.set_xlabel(xlab)
    if ylab == None:
        pass
    else:
        ax.set_ylabel(ylab)
        
    # Check if grid lines are desired
    if grid_on:
        ax.grid()
       
    # Set the colorbar if desired
    if cb_flag:
        cb = plt.colorbar(mappable=p,orientation=cb_orient)
        if cb_lab == None:
            pass
        else:
            cb.set_label(cb_lab)
            
    return p
#**===============================================================
def sweep_to_Cart(Var,range,rot,tilt,ax=None,data_proj='fore',title=None,
               vmin=-24.,vmax=64.,cmap=None,mask_outside=True,grid_on=True,
               xlims=None,ylims=None,xlab=None,ylab=None,
               cb_flag=True,cb_orient='horizontal',cb_lab=None):
    """Project native (polar) coordinate radar sweep data onto a flat Cartesian coordinate grid.
       See Lee et al. (1994) Journal of Atmospheric and Oceanic Technology for methodology and definitions.
 INPUT::
  Var             = Data values to plot
  range           = Range along ray
  rot             = Rotation angle with respect to instrument [degrees]
  tilt            = Tilt angle with respect to platform [degrees]

       THESE NEXT VALUES CUSTOMIZE THE PLOT
  data_proj       = Which direction the data is collected [ 'fore' or 'aft' ]
                    Needed to display the data in forward-looking convention
  vmin            = Minimum value to display
  vmax            = Maximum value to display
  cmap            = Matplotlib colormap name
  mask_outside    = Boolean flag to mask values outside vmin/vmax
  title           = Title to label plot with, if none given it is omitted
  xlims           = X-axis limits (min,max)
  ylims           = Y-axis limits (min,max)
  cb_flag         = True turns on colorbar, False no colorbar
  cb_lab          = Colorbar label, None will use a default label generated from the
                     field information.
  cb_orient       = Orientation of colorbar (vertical or horizontal)
 OUTPUT::
  fig             = Figure to create plot to
  ax              = Polar axis created
  p               = Output plot
 USAGE::
  p, fig, ax = DPJgrid_Horiz(fig,Z,W,Lon,Lat,Ht,TLon,TLat,TAlt,[Zcoord],[Ucoord],[WindVec],[Track],
                     [cminmax],[clevs],[vmin],[vmax],[dlat],[dlon],[proj],
                     [cmap],[title],[pName],[pType],[figsize])
 NOTES::
 Plotting convention is to project the data onto a 2D surface looking from the back of
   aircraft forward.  The polar data has no direct measure of direction and therefore for 
   the P-3, determination of 'fore' or 'aft' may be accomplished by looking at tilt angle.
   Positive is fore, negative is aft.
    """
# HISTORY::
#  15 Apr 2014 - Nick Guy NOAA/NSSL/WRDD, NRC 
#--------------------------------------------------------
    # Check the data projection
    if data_proj == 'fore':
        Fact = 1.
    elif data_proj == 'aft':
        Fact = -1.
    else:
        print "Need the data collection direction, either 'fore' or 'aft', assuming 'fore'"
        Fact = 1.
        
    range = np.array(range) # Make sure that the data is numpy array
 
    values = np.array(Var) # Make sure that the data is numpy array
    values = values.reshape(len(rot), len(range)) # Resize it to work
    
    # mask the data where outside the limits
    if mask_outside:
        Var = np.ma.masked_outside(Var, vmin, vmax)

    # Create 2D variables to plot contour against
    r, Rot = np.meshgrid(range, np.radians(rot))
    r2, Tilt = np.meshgrid(range,np.radians(tilt))
    
    # Convert r, Rot to Cartesian (x,z)
    X = Fact * r * np.sin(Rot) * np.sin(Tilt)
    Z = r * np.cos(Rot) * np.cos(Tilt)
    
    # Set axis limits if passed
    if xlims == None:
        xlims = (-1.05 * range.max(), 1.05 * range.max())
    else:
        xlims = xlims
    if ylims == None:
        ylims = (-10.,30.)
    else:
        ylims = ylims
    
    # Plot the data
    p = plot_sweep_grid(X,Z,Var,ax=ax,title=title,
               vmin=vmin,vmax=vmax,cmap=cmap,grid_on=grid_on,
               xlims=xlims,ylims=ylims,xlab=xlab,ylab=ylab,
               cb_flag=cb_flag,cb_orient=cb_orient,cb_lab=cb_lab)
            
    return p
#**===============================================================
def sweep_aircraft_relative(Var,range,tilt,rot,ax=None,title=None,
               vmin=-24.,vmax=64.,cmap=None,mask_outside=True,grid_on=True,
               xlims=None,ylims=None,xlab=None,ylab=None,
               cb_flag=True,cb_orient='horizontal',cb_lab=None):
    """Project native (polar) coordinate radar sweep data onto aircraft-relative Cartesian coordinate grid.
       See Lee et al. (1994) Journal of Atmospheric and Oceanic Technology for methodology and definitions.
 INPUT::
  Var             = Data values to plot
  range           = Range along ray
  tilt            = Radar ray Tilt angle with respect to platform [degrees]
  rot             = Rotation angle with respect to instrument [degrees]

       THESE NEXT VALUES CUSTOMIZE THE PLOT
  vmin            = Minimum value to display
  vmax            = Maximum value to display
  cmap            = Matplotlib colormap name
  mask_outside    = Boolean flag to mask values outside vmin/vmax
  title           = Title to label plot with, if none given it is omitted
  xlims           = X-axis limits (min,max)
  ylims           = Y-axis limits (min,max)
  cb_flag         = True turns on colorbar, False no colorbar
  cb_lab          = Colorbar label, None will use a default label generated from the
                     field information.
  cb_orient       = Orientation of colorbar (vertical or horizontal)
 OUTPUT::
  fig             = Figure to create plot to
  ax              = Polar axis created
  p               = Output plot
 USAGE::
  p = sweep_aircraft(Var,range,tilt,rot,[**args])
 NOTES::
  Plotting convention is to project the data onto a 2D surface looking from the back of
   aircraft forward.  X,Y,Z coordinates are a direct projection from polar to Cartesian
   coordinates.
  This mapping does NOT take into account corrections for roll, pitch, or drift of the
   aircraft.
    """
# HISTORY::
#  15 Apr 2014 - Nick Guy NOAA/NSSL/WRDD, NRC 
#--------------------------------------------------------
    range = np.array(range) # Make sure that the data is numpy array
 
    values = np.array(Var) # Make sure that the data is numpy array
    values = values.reshape(len(rot), len(range)) # Resize it to work
    
    # mask the data where outside the limits
    if mask_outside:
        Var = np.ma.masked_outside(Var, vmin, vmax)

    # Create 2D variables to plot contour against
    r, Rot = np.meshgrid(range, np.radians(rot))
    r2, Tilt = np.meshgrid(range,np.radians(tilt))
    
    # Convert polar (r,Rot,Tilt) to Cartesian (x,y,z)
    X = r * np.cos(Tilt) * np.sin(Rot)
    Y = r * np.sin(Tilt)
    Z = r * np.cos(Rot) * np.cos(Tilt)
    
    # Set axis limits if passed
    if xlims == None:
        xlims = (-1.05 * range.max(), 1.05 * range.max())
    else:
        xlims = xlims
    if ylims == None:
        ylims = (-10.,30.)
    else:
        ylims = ylims
    
    # Plot the data
    p = plot_sweep_grid(X,Z,Var,ax=ax,title=title,
               vmin=vmin,vmax=vmax,cmap=cmap,grid_on=grid_on,
               xlims=xlims,ylims=ylims,xlab=xlab,ylab=ylab,
               cb_flag=cb_flag,cb_orient=cb_orient,cb_lab=cb_lab)

    return p
#**===============================================================
def sweep_track_relative(Var,range,tilt,rot,roll,drift,pitch,ax=None,title=None,
               vmin=-24.,vmax=64.,cmap=None,mask_outside=True,grid_on=True,
               xlims=None,ylims=None,xlab=None,ylab=None,
               cb_flag=True,cb_orient='horizontal',cb_lab=None):
    """Project native (polar) coordinate radar sweep data onto track-relative Cartesian coordinate grid.
       See Lee et al. (1994) Journal of Atmospheric and Oceanic Technology for methodology and definitions..
 INPUT::
  Var             = Data values to plot
  range           = Range along ray
  tilt            = Radar ray Tilt angle with respect to platform [degrees]
  rot             = Rotation angle with respect to instrument [degrees]
  roll            = Aircraft roll angle [degrees], right-wing down positive
  drift           = Drift angle [degrees], between heading and track
  pitch           = Pitch angle [degrees], nose up positive

       THESE NEXT VALUES CUSTOMIZE THE PLOT
  vmin            = Minimum value to display
  vmax            = Maximum value to display
  cmap            = Matplotlib colormap name
  mask_outside    = Boolean flag to mask values outside vmin/vmax
  title           = Title to label plot with, if none given it is omitted
  xlims           = X-axis limits (min,max)
  ylims           = Y-axis limits (min,max)
  cb_flag         = True turns on colorbar, False no colorbar
  cb_lab          = Colorbar label, None will use a default label generated from the
                     field information.
  cb_orient       = Orientation of colorbar (vertical or horizontal)
 OUTPUT::
  fig             = Figure to create plot to
  ax              = Polar axis created
  p               = Output plot
 USAGE::
  p = sweep_track_relative(Var,range,tilt,rot,roll,drift,pitch,[**args])
 NOTES::
  Plotting convention is to project the data onto a 2D surface looking from the back of
   aircraft forward.  X,Y,Z coordinates are a rotation of the data about the aircraft
   track, following a direct projection from polar to Cartesian coordinates.
  This mapping corrects for roll, pitch, and drift of the aircraft.
  This is considered a leveled, heading-relative coordinate system.
    """
# HISTORY::
#  15 Apr 2014 - Nick Guy NOAA/NSSL/WRDD, NRC 
#--------------------------------------------------------
    range = np.array(range) # Make sure that the data is numpy array
 
    values = np.array(Var) # Make sure that the data is numpy array
    values = values.reshape(len(rot), len(range)) # Resize it to work
    
    # mask the data where outside the limits
    if mask_outside:
        Var = np.ma.masked_outside(Var, vmin, vmax)

    # Create 2D variables to plot contour against
    r, Rot = np.meshgrid(range, np.radians(rot))
    r2, Tilt = np.meshgrid(range,np.radians(tilt))
    r3, Roll = np.meshgrid(range,np.radians(roll))
    r4, Pitch = np.meshgrid(range,np.radians(pitch))
    r5, Drift = np.meshgrid(range,np.radians(drift))
    del r2,r3,r4,r5
    
    # Convert r, Rot to Cartesian (x,y)
    X = r * (np.cos(Rot + Roll) * np.sin(Drift) * np.cos(Tilt) * np.sin(Pitch) +
        np.cos(Drift) * np.sin(Rot + Roll) * np.cos(Tilt) - 
        np.sin(Drift) * np.cos(Pitch) * np.sin(Tilt))
    Y = r * (-1.* np.cos(Rot + Roll) * np.cos(Drift) * np.cos(Tilt) * np.sin(Pitch) +
        np.sin(Drift) * np.sin(Rot + Roll) * np.cos(Tilt) + 
        np.cos(Drift) * np.cos(Pitch) * np.sin(Tilt))
    Z = r * np.cos(Pitch) * np.cos(Tilt) * np.cos(Rot + Roll) + np.sin(Pitch) * np.sin(Tilt)
    
    # Set axis limits if passed
    if xlims == None:
        xlims = (-1.05 * range.max(), 1.05 * range.max())
    else:
        xlims = xlims
    if ylims == None:
        ylims = (-10.,30.)
    else:
        ylims = ylims
    
    # Plot the data
    p = plot_sweep_grid(X,Z,Var,ax=ax,title=title,
               vmin=vmin,vmax=vmax,cmap=cmap,grid_on=grid_on,
               xlims=xlims,ylims=ylims,xlab=xlab,ylab=ylab,
               cb_flag=cb_flag,cb_orient=cb_orient,cb_lab=cb_lab)
            
    return p
#**===============================================================
def sweep_earth_relative(Var,range,tilt,rot,roll,heading,pitch,ax=None,title=None,
               vmin=-24.,vmax=64.,cmap=None,mask_outside=True,grid_on=True,
               xlims=None,ylims=None,xlab=None,ylab=None,
               cb_flag=True,cb_orient='horizontal',cb_lab=None):
    """Project native (polar) coordinate radar sweep data onto earth-relative Cartesian coordinate grid.
       See Lee et al. (1994) Journal of Atmospheric and Oceanic Technology for methodology and definitions.
 INPUT::
  Var             = Data values to plot
  range           = Range along ray
  tilt            = Radar ray Tilt angle with respect to platform [degrees]
  rot             = Rotation angle with respect to instrument [degrees]
  roll            = Aircraft roll angle [degrees], right-wing down positive
  heading         = Heading angle [degrees], clockwise from North
  pitch           = Pitch angle [degrees], nose up positive

       THESE NEXT VALUES CUSTOMIZE THE PLOT
  vmin            = Minimum value to display
  vmax            = Maximum value to display
  cmap            = Matplotlib colormap name
  mask_outside    = Boolean flag to mask values outside vmin/vmax
  title           = Title to label plot with, if none given it is omitted
  xlims           = X-axis limits (min,max)
  ylims           = Y-axis limits (min,max)
  cb_flag         = True turns on colorbar, False no colorbar
  cb_lab          = Colorbar label, None will use a default label generated from the
                     field information.
  cb_orient       = Orientation of colorbar (vertical or horizontal)
 OUTPUT::
  fig             = Figure to create plot to
  ax              = Polar axis created
  p               = Output plot
 USAGE::
  p = sweep_earth_relative(Var,range,tilt,rot,roll,heading,pitch,[**args])
 NOTES::
  Plotting convention is to project the data onto a 2D surface looking from the back of
   aircraft forward.  X,Y,Z coordinates are a rotation of the data about an 
   earth-relative azimuth, following a direct projection from polar to 
   Cartesian coordinates.
  This mapping corrects for roll, pitch, and drift of the aircraft.
  This is considered a leveled, heading-relative coordinate system.
    """
# HISTORY::
#  15 Apr 2014 - Nick Guy NOAA/NSSL/WRDD, NRC 
#--------------------------------------------------------
    range = np.array(range) # Make sure that the data is numpy array
 
    values = np.array(Var) # Make sure that the data is numpy array
    values = values.reshape(len(rot), len(range)) # Resize it to work
    
    # mask the data where outside the limits
    if mask_outside:
        Var = np.ma.masked_outside(Var, vmin, vmax)

    # Create 2D variables to plot contour against
    r, Rot = np.meshgrid(range, np.radians(rot))
    r2, Tilt = np.meshgrid(range,np.radians(tilt))
    r3, Roll = np.meshgrid(range,np.radians(roll))
    r4, Pitch = np.meshgrid(range,np.radians(pitch))
    r5, Heading = np.meshgrid(range,np.radians(heading))
    del r2,r3,r4,r5
    
    # Convert r, Rot to Cartesian (x,y)
    X = r * (-1.* np.cos(Rot + Roll) * np.sin(Heading) * np.cos(Tilt) * np.sin(Pitch) +
        np.cos(Heading) * np.sin(Rot + Roll) * np.cos(Tilt) + 
        np.sin(Heading) * np.cos(Pitch) * np.sin(Tilt))
    Y = r * (-1.* np.cos(Rot + Roll) * np.cos(Heading) * np.cos(Tilt) * np.sin(Pitch) -
        np.sin(Heading) * np.sin(Rot + Roll) * np.cos(Tilt) + 
        np.cos(Heading) * np.cos(Pitch) * np.sin(Tilt))
    Z = r * np.cos(Pitch) * np.cos(Tilt) * np.cos(Rot + Roll) + np.sin(Pitch) * np.sin(Tilt)
    
    # Set axis limits if passed
    if xlims == None:
        xlims = (-1.05 * range.max(), 1.05 * range.max())
    else:
        xlims = xlims
    if ylims == None:
        ylims = (-10.,30.)
    else:
        ylims = ylims
    
    # Plot the data
    p = plot_sweep_grid(X,Z,Var,ax=ax,title=title,
               vmin=vmin,vmax=vmax,cmap=cmap,grid_on=grid_on,
               xlims=xlims,ylims=ylims,xlab=xlab,ylab=ylab,
               cb_flag=cb_flag,cb_orient=cb_orient,cb_lab=cb_lab)
    return p
#**===============================================================