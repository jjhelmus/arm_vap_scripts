#!/usr/bin/env python

import matplotlib.pyplot as plt
import netCDF4
import pyart

# MMCG figure
dataset = netCDF4.Dataset('sgpcsaprmmcgi7.c0.20110520.110100.nc')
refl = dataset.variables['reflectivity_horizontal']
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(refl[0,4], origin='lower')
fig.savefig('mapped_figure.png')

# Test
dataset = netCDF4.Dataset('foo.dir/sgpcsaprmmcgI7.c0.20110520.110100.nc')
refl = dataset.variables['reflectivity_horizontal']
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(refl[0,4], origin='lower')
fig.savefig('exp_figure.png')


# Radial coords
"""
radar = pyart.io.read_netcdf('sgpcsaprsurcmacI7.c0.20110520.110100.nc')
display = pyart.graph.RadarDisplay(radar)
fig = plt.figure()
ax = fig.add_subplot(111)
display.plot_ppi('reflectivity_horizontal', 0, vmin=-16, vmax=48)
fig.savefig('radial_figure.png')
"""
