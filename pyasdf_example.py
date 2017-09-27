import wfdiff.specfem_helper as sh
import pyasdf
import obspy

iex = 1
# 1: Example of creating an asdf file
# 2: Example of reading an asdf file

# Convert specfem file to asdf
if iex == 1:
    # Need station and event files
    stations_file = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/STATIONS_FILTERED'
    event_file = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/CMTSOLUTION'

    # Add one low res file
    folder = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/'
    asdf_filename = '/home/vipul/REPOSITORIES/ADJOINT_TOMO/wfdiff_uaf/nenana_gll5.h5' 
    wf_tag = 'gll5'
    sh.convert_to_asdf(asdf_filename, folder, stations_file, event_file, wf_tag)

    # Also add a high res file for comparison
    folder = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL7/'
    asdf_filename = '/home/vipul/REPOSITORIES/ADJOINT_TOMO/wfdiff_uaf/nenana_gll7.h5' 
    wf_tag = 'gll7'
    sh.convert_to_asdf(asdf_filename, folder, stations_file, event_file, wf_tag)

# Read asdf as stream
if iex == 2:
    asdf_filename = 'nenana.h5'
    wf_tag = 'gll5'
    ds = pyasdf.ASDFDataSet(asdf_filename)
    st_low = sh.get_stream_from_asdf(ds,wf_tag)
    print(st_low.__str__(extended=True))

    # Read high res files
    wf_tag = 'gll7'
    st_high = sh.get_stream_from_asdf(ds,wf_tag)
    print(st_high.__str__(extended=True))
