# python specfem_util_test.py

import obspy
import wfdiff.specfem_helper as sh

import os

cmtsolution_file = 'test_data/CMTSOLUTION'
station_file = 'test_data/STATIONS'
waveforms = 'test_data/NGLL5/*semd'
# Also run this
# To create RTZ for NGLL7 case
waveforms = 'test_data/NGLL7/*semd'


#-----------------------------------------------------------------

# Read CMTSOLUTION file
ev = sh.read_specfem_cmtsolution_file(cmtsolution_file)
eid = sh.otime2eid(ev.origins[0].time)

# Read station file
stn = sh.read_specfem_stations_file(station_file)

# add event info
st = sh.add_event_info(ev, stn, waveforms)

# save ENZ as sac
ENZ_dir = eid +'_ENZ'
sh.save_as_sac(st, ENZ_dir)

# Rotate to RTZ
print('---> Rotating to RTZ')
st.rotate('NE->RT')

# save RTZ as sac
RTZ_dir = eid +'_RTZ'
sh.save_as_sac(st, RTZ_dir)
