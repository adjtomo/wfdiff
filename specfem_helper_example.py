# python specfem_util_test.py

import obspy
import wfdiff.specfem_helper as sh

import os

event_file = 'test_data/CMTSOLUTION'
stations_file = 'test_data/STATIONS'
specfem_wf = 'test_data/NGLL5/*semd'
# Also run this
# To create RTZ for NGLL7 case
#specfem_wf = 'test_data/NGLL7/*semd'

# Nenana example
stations_file = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/STATIONS_FILTERED'
event_file = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/CMTSOLUTION'
specfem_wf = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/*semd'

#-----------------------------------------------------------------

# Read CMTSOLUTION file
cat = obspy.read_events(event_file)
event = cat[0]
eid = event.origins[0].time.strftime("%Y%m%d%H%M%S%f")[:-3]

# Read specfem station file as pandas dataframe
stn = sh.read_specfem_stations_file(stations_file)

# Read specfem files
st = sh.read_specfem_files(specfem_wf)

# Add backazimuth. Needed for rotation
st = sh.add_event_station_info(st, event, stn)

# save ENZ as sac
# NOTE: trace.stats.back_azimuth will not be saved to sac
# Put in sac header if you want to have it saved
# And repopulate the trace.stat when you reload the sac files
ENZ_dir = eid +'_ENZ'
sh.save_as_sac(st, ENZ_dir)

#-------------------------------------------------
# Convert sac to asdf
from pyasdf.scripts import sac2asdf

asdf_filename = eid + '.h5'
sac_dir = ENZ_dir
tag = 'enz'
# NOTE: These will not have any backazimuth info
# Best way is to skip the SAC step and directly save traces
# as asdf. This will save back_azimuth and other trace.stats attributes
sac2asdf.add_to_adsf_file('test_asdf.h5', sac_dir, tag)

#-------------------------------------------------
# Rotate to RTZ
print('---> Rotating to RTZ')
st.rotate('NE->RT')

# save RTZ as sac
RTZ_dir = eid +'_RTZ'
sh.save_as_sac(st, RTZ_dir)

