#! /usr/bin/env python
# sample script for processing phase using LP methods
# usage : process_and_save_xsar_params.py filename outdir params

import sys
import os
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir)

import matplotlib as mpl
mpl.use('Agg')
import copy
import netCDF4
import pyart
from pyart.io import radar, nc_utils
from pyart.correct import phase_proc, attenuation, dealias
import getpass
import socket
import numpy as np


def dt_to_dict(dt, **kwargs):
        pref = kwargs.get('pref', '')
        return dict([(pref+key, getattr(dt, key)) for key in
                    ['year', 'month', 'day', 'hour', 'minute', 'second']])


def multiconvert(item, typ):
    if typ == 'string':
        ret = item.strip()
    if typ == 'bool':
        if item.strip().lower() == 'true':
            ret = True
        else:
            ret = False
    if typ == 'float':
        ret = float(item)
    if typ == 'int':
        ret = int(item)
    return ret


def parse_prefs(filename):
    bf = open(filename, 'r')
    text = bf.readlines()
    bf.close()
    params = {}
    metadata = {}
    for item in text:
        details = item.split(' ')
        data = item.split(details[1])
        if 'metadata' in details[0]:
            metadata.update({details[0].split('-')[0]: multiconvert(data[1],
                            details[1])})
        else:
            params.update({details[0]: multiconvert(data[1], details[1])})
    return params, metadata

if __name__ == "__main__":
    process_version="1.0E"
    dod_version="1.0"
    # parse command line arguments
    filename = sys.argv[1]
    outdir = sys.argv[2]
    params, metadata = parse_prefs(sys.argv[3])
    is_dir = params['sounding_dir']

    # read in radar data
    myradar = pyart.io.read_rsl(filename)
    #myradar = radar.Radar(py4dd.RSL_anyformat_to_radar(filename))
    #myradar.metadata.update(metadata)

    # dealias
    target = netCDF4.num2date(myradar.time['data'][0], myradar.time['units'])
    fname = ('sgpinterpolatedsondeC1.c1.%(year)04d%(month)02d%(day)02d.000000.cdf' % nc_utils.dt_to_dict(target))
    print fname
    interp_sonde = netCDF4.Dataset(is_dir + fname)
    myheight, myspeed, mydirection = dealias.find_time_in_interp_sonde(
        interp_sonde, target)
    my_new_field = pyart.correct.dealias_fourdd(myradar, myheight * 1000.0, myspeed,
                                 mydirection, target)
    myradar.fields.update({'corrected_mean_doppler_velocity': my_new_field})
    interp_sonde.close()
    print("done")
    # process Phase
    gates = myradar.range['data'][1] - myradar.range['data'][0]
    rge = 20.0 * 1000.0
    ng = rge / gates
    # append a datetime object
    mydatetime = netCDF4.num2date(myradar.time['data'][0],
                                  myradar.time['units'],
                                  calendar=myradar.time['calendar'])
    mydict = dt_to_dict(mydatetime)
    mydict.update({'scanmode': {'ppi': 'sur', 'rhi': 'rhi'}
                  [myradar.sweep_mode[0]], 'fac': metadata['facility'],
                  'name':metadata['instrument_name']})
    ofilename = outdir + '%(name)s%(scanmode)scmac%(fac)s.c0.%(year)04d%(month)02d%(day)02d.%(hour)02d%(minute)02d%(second)02d.nc' % mydict
    reproc_phase, sob_kdp = phase_proc.phase_proc_lp(
        myradar, params['reflectivity_offset'], sys_phase=params['sys_phase'],
        overide_sys_phase=params['overide_sys_phase'], debug=True, nowrap=ng)
    myradar.fields.update({'recalculated_diff_phase': sob_kdp,
                         'proc_dp_phase_shift': reproc_phase})
    #attenuation correction
    spec_at, cor_z = attenuation.calculate_attenuation(
       myradar, params['reflectivity_offset'], debug=True, a_coef=0.17)
    myradar.fields.update({'specific_attenuation': spec_at})
    myradar.fields.update({'corrected_reflectivity_horizontal': cor_z})
    R = 51.3 * (myradar.fields['specific_attenuation']['data']) ** 0.81
    rainrate = copy.deepcopy(myradar.fields['diff_phase'])
    rainrate['data'] = R
    rainrate['valid_min'] = 0.0
    rainrate['valid_max'] = 400.0
    rainrate['standard_name'] = 'rainfall_rate'
    rainrate['long_name'] = 'rainfall_rate'
    rainrate['least_significant_digit'] = 1
    rainrate['units'] = 'mm/hr'
    myradar.fields.update({'rain_rate_A': rainrate})
    mask=myradar.fields['corrected_reflectivity_horizontal']['data'].mask
    myradar.fields['rain_rate_A']['data'][np.where(mask)]=0.0
    myradar.fields['rain_rate_A'].update({'comment':'Rain rate calculated from specific_attenuation, R=51.3*specific_attenuation**0.81, note R=0.0 where norm coherent power < 0.4 or rhohv < 0.8'})
    line=""
    for arg in sys.argv:
        line = line + arg + " "
    myradar.metadata.update({'command_line':line, 'process_version':process_version, 'dod_version':dod_version})
    pyart.io.write_netcdf(ofilename, myradar, format='NETCDF4')
    #netcdf_obj = netCDF4.Dataset(ofilename, 'w', format='NETCDF4')
    #nc_utils.write_radar4(netcdf_obj, myradar)
    #netcdf_obj.close()
