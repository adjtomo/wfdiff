import wfdiff.specfem_helper as sh
import pyasdf
import obspy

iex = 0
# 0: Default example
# 1: Example of creating an asdf file
# 2: Example of reading an asdf file
# 3: Read asdf file rotate to RTZ

if iex == 0:
    # Need station and event files
    stations_file = "./test_data/STATIONS"
    event_file = './test_data/CMTSOLUTION'

    # Convert low res file to asdf
    specfem_folder = "./test_data/NGLL5"
    asdf_filename = "./test_data/nenana_gll5.h5"
    wf_tag = 'gll5'
    sh.specfem_to_asdf(asdf_filename, specfem_folder, stations_file, event_file, wf_tag)

    # Convert low res file to asdf
    specfem_folder = "./test_data/NGLL7"
    asdf_filename = "./test_data/nenana_gll7.h5"
    wf_tag = 'gll7'
    sh.specfem_to_asdf(asdf_filename, specfem_folder, stations_file, event_file, wf_tag)


# Convert specfem file to asdf
if iex == 1:
    # Need station and event files
    stations_file = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/STATIONS_FILTERED'
    event_file = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/CMTSOLUTION'

    # Add one low res file
    specfem_folder = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/'
    asdf_filename = '/home/vipul/REPOSITORIES/ADJOINT_TOMO/wfdiff_uaf/nenana_gll5.h5' 
    wf_tag = 'gll5'
    sh.specfem_to_asdf(asdf_filename, specfem_folder, stations_file, event_file, wf_tag)

    # Also add a high res file for comparison
    #specfem_folder = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL7/'
    #asdf_filename = '/home/vipul/REPOSITORIES/ADJOINT_TOMO/wfdiff_uaf/nenana_gll7.h5' 
    #wf_tag = 'gll7'
    #sh.specfem_to_asdf(asdf_filename, specfem_folder, stations_file, event_file, wf_tag)

# Read asdf as stream
if iex == 2:
    asdf_filename = 'nenana_gll5.h5'
    wf_tag = 'gll5'
    ds = pyasdf.ASDFDataSet(asdf_filename)
    st_low = sh.get_stream_from_asdf(ds,wf_tag)
    print(st_low.__str__(extended=True))

    # Read high res files
    asdf_filename = 'nenana_gll7.h5'
    wf_tag = 'gll7'
    ds = pyasdf.ASDFDataSet(asdf_filename)
    st_high = sh.get_stream_from_asdf(ds,wf_tag)
    print(st_high.__str__(extended=True))

# Read stream from asdf file and rotate
if iex == 3:
    irotate = True

    asdf_filename = 'nenana_gll5.h5'
    wf_tag = 'gll5'
    print('Reading asdf dataset')
    ds = pyasdf.ASDFDataSet(asdf_filename)

    print('Extracting streams, event and station info')
    st_ENZ = sh.get_stream_from_asdf(ds,wf_tag)
    event = ds.events[0]
    stations = sh.get_stations_from_asdf(ds)
    # Now add back_azimuth 
    if irotate is True:
        st = sh.add_event_station_info(st_ENZ, event, stations)
        print('Rotating')
        st.rotate('NE->RT')
        ds.add_waveforms(st, tag='rtz', event_id=event)
    print(st)
