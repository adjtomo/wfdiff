# python specfem_util_test.py

import obspy
import wfdiff.specfem_helper as sh
import shutil
import os

event_file = '../CMTSOLUTION'
stations_file = '../STATIONS'
specfem_wf = './*semd'
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
st.resample(50)
# Rotate to RTZ
print('---> Rotating to RTZ')
st.rotate('NE->RT')

# save RTZ as sac
RTZ_dir = eid +'_RTZ'
sh.save_as_sac(st, RTZ_dir)

#files = os.listdir()
#for file in files:
#    a = file[:-4]
#    b = a[-1]
#    c = 'Mtt'
#    print(a[:-3]+'.'+b+'.'+c+'.sac')
#    shutil.move(file, a[:-3]+'.'+b+'.'+c+'.sac')
