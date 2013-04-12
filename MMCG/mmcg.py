#!/usr/bin/env python
#mapping rainmaps from CSAPR CMAC data
import matplotlib.pyplot as plt

import sys

import netCDF4
import numpy as np

import pyart
from pyart.io import grid, nc_utils

# requirements which should be removed????
from pyart.io.common import dms_to_d
from pyart.graph.common import corner_to_point
from pyart.io.nc_utils import dt_to_dict


def my_qrf(xg, yg, zg):
    hz_dis = np.sqrt(xg ** 2 + yg ** 2)
    hyp = np.sqrt(hz_dis ** 2 + zg ** 2) + 1.0
    theta = np.arccos(hz_dis / hyp) * 180.0 / np.pi
    qrf = (zg / 20.0) + np.sqrt(yg ** 2 + xg ** 2) * \
        np.tan(1.0 * 1.5 * np.pi / 180.0) + 500.0 + \
        zg * 1000.0 * theta ** 2 / (45.0 ** 2) / 17000.
    return qrf


if __name__ == "__main__":
    process_version = "1.0E"
    vapname = "mmcg"
    site = "sgp"
    facility = "I7"
    level = "c0"
    ncfile = sys.argv[1]
    odir = sys.argv[2]
    cf_alt = 320.0 / 1000.0
    cf_lat = dms_to_d([36.0, 36.0, 18.35])
    cf_lon = -1.0*dms_to_d([97.0, 29.0, 10.69])
    #mync=netCDF4.Dataset(ncfile)
    #myradar=radar.Radar(mync)
    myradar = pyart.io.read_netcdf(ncfile)
    print myradar.time['data'][0]
    cp = corner_to_point([myradar.location['latitude']['data'],
                          myradar.location['longitude']['data']],
                         [cf_lat, cf_lon])
    #mync.close()
    mygrids = grid.pyGrid(
        (myradar,), nxyz=(241, 241, 35),
        xyzr=((-120000 - cp[0], 120000 - cp[0]),
             (-120000 - cp[1], 120000 - cp[1]),
             (0, 17000)),
        params=['reflectivity_horizontal'],
        toa=20000, origin=[cf_lat, cf_lon, cf_alt], qrf=my_qrf)

    # plot for testing...
    refl = mygrids.fields['reflectivity_horizontal']['data']
    refl = np.ma.masked_equal(refl, -9999.0)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(refl[4], origin='lower')
    fig.savefig('foo_figure.png')

    mydatetime = netCDF4.num2date(myradar.time['data'][0],
                                  myradar.time['units'],
                                  calendar=myradar.time['calendar'])
    time_str = mydatetime.strftime('%Y%m%d.%H%M%S')

    desc = """Raw and corrected moments from the SGP CSAPR
            scanning radar mapped to a cartesian grid.
            The origin (x=y=z=0 is at the SGP Central facility.
            Point were mapped using a Barnes type filter with
            a radius of influence (ROI) that varies as:
            ROI= (z / 20.0) + np.sqrt(y ** 2 + x ** 2) * \
            np.tan(1.0 * 1.5 * np.pi / 180.0) + 500.0  + \
            zg*1000.0*antenna_elevation**2/(45.0**2)/17000."""

    mygrids.metadata.update({'description': desc})
    mygrids.metadata.update({'institution': 'ARM Climate research facility'})
    mygrids.metadata.update(
        {'contact':
            'Scott Collis, Argonne National Laboratory, scollis@anl.gov'})
    mygrids.metadata['input_datastreams_num'] = "1"
    mygrids.metadata.update({'input_datastreams': ncfile.split('/')[-1]})
    mygrids.metadata.pop('Conventions')
    mygrids.metadata.pop('source')
    t = (myradar.metadata['instrument_name'], vapname,
         myradar.metadata['facility'], level, time_str)
    myfilename = odir + "%s%s%s.%s.%s.nc" % t
    ncfobj = netCDF4.Dataset(myfilename, 'w', format='NETCDF4')
    nc_utils.save_pyGrid(ncfobj, mygrids, zlib=True)
    ncfobj.close()
    print mygrids.axes['time_start']['data']
