# To run this script
# mpirun -n 8 python run_wfdiff_test.py

from wfdiff.wfdiff import WFDiff

# A non-existant directory where the output will be stored.
#output_directory = "output"

# Choose which misfits you want and the threshold values.
THRESHOLDS = {
    "rms": 0.25,
    #"l1_norm": 0.3,
    # A couple more misfits are available.
    # These currently take a long time to calculate.
    #"cross_correlation": 0.9,
    #"phase_misfit": 0.1,
    #"envelope_misfit": 0.1
    }

#--------------------------------------------------------------------------------
# Run example
event_file = None
stations_file = None

low_res_seismos = None    # low resolution files
high_res_seismos = None   # high resolution files
wf_format = None          # these could be sac, specfem, asdf or any other obspy suppported format
trace_tags = None         # tags for low_res_seismos and high_res_seismos


iex = 0
# 0: Default example

elif iex == 0:
    # Use asdf file instead (Nenana example)
    low_res_seismos="./test_data/nenana_gll5.h5"
    high_res_seismos="./test_data/nenana_gll7.h5"
    wf_format = 'asdf'
    trace_tags = ['gll5','gll7']  # for waveform debug plots
    asdf_tags = ['gll5','gll7']
    output_directory = "output_test"
    output_format = 'pdf' # or ps, eps etc.

if iex == 1:
    # test data (using ENZ components)
    low_res_seismos="./test_data/NGLL5/*semd"
    high_res_seismos="./test_data/NGLL7/*semd"
    stations_file="./test_data/STATIONS"
    event_file='./test_data/CMTSOLUTION'
    trace_tags = ['NGLL5','NGLL7']
    wf_format = 'specfem'
    output_directory = "output_test"
    output_format = 'pdf' # or ps, eps etc.

elif iex == 2: 
    # Test using sac data (using ENZ components) - Same as the Example 1
    # Run specfem_util_test.py first to create rotate SAC file for NGLL5 and NGLL7 case
    # test (using ENZ components) - Same as the Example 1
    low_res_seismos="./20140831030657110_NGLL5_ENZ/*sac"
    high_res_seismos="./20140831030657110_NGLL7_ENZ/*sac"
    trace_tags = ['NGLL5-ENZ','NGLL7-ENZ']
    wf_format = 'sac'
    output_directory = "output_test"
    output_format = 'pdf' # or ps, eps etc.

elif iex == 3:
    # Test using sac data (using RTZ components)
    # Run specfem_util_test.py first to create rotate SAC file for NGLL5 and NGLL7 case
    low_res_seismos="./20140831030657110_NGLL5_RTZ/*sac"
    high_res_seismos="./20140831030657110_NGLL7_RTZ/*sac"
    stations_file="./test_data/STATIONS"
    trace_tags = ['NGLL5-RTZ','NGLL7-RTZ']
    wf_format = 'sac'
    output_directory = "output_test"
    output_format = 'pdf' # or ps, eps etc.

elif iex == 4:
    # Use asdf file instead (Nenana example)
    # See pyasdf_example.py for creating an asdf file from specfem output files 
    # 
    #low_res_seismos="/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/*.semd"
    #high_res_seismos="/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL7/*.semd"
    #stations_file = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/STATIONS_FILTERED'
    #event_file = '/store/homeglut/carltape/OUTPUT_SPECFEM3D/alaska/nenana/OUTPUT_FILES_nenana_NGLL5/CMTSOLUTION'
    #wf_format = 'specfem'
    low_res_seismos="./nenana_gll5.h5"
    high_res_seismos="./nenana_gll7.h5"
    wf_format = 'asdf'
    trace_tags = ['gll5','gll7']  # for waveform debug plots
    asdf_tags = ['gll5','gll7']
    output_directory = "output_test"
    output_format = 'pdf' # or ps, eps etc.
#--------------------------------------------------------------------------------

# Configuration.
c = WFDiff(
    low_res_seismos=low_res_seismos, high_res_seismos=high_res_seismos, 
    stations_file=stations_file, event_file=event_file,
    # Specify the units of the data and the units the analysis should take
    # place in.
    data_units="displacement", desired_analysis_units="displacement",
    # Periods to test.
    t_min=.1, t_max=1, dt=.1,
    # Data window to take into account in seconds since the first sample.
    starttime=20, endtime=100,
    rotate_RTZ = True,
    # Set to 'True' if specfem filename are in NET.STA.CHA.* format
    # Set to 'False' if specfem filename are in STA.NET.CHA.* format
    new_specfem_name=True,
    # Tags printed on the plots
    # Note: for the asdf files these are same as the waveform dataset tag
    trace_tags=trace_tags,
    asdf_tags=asdf_tags,
    # Set to 'specfem' if waveform ASCII files are used. 
    # Set to 'asdf' if asdf waveform file
    # All other fileformat should otherwise work just fine.
    wf_format=wf_format)


# Perform the frequency dependent misfit measurements. Results will be stored
# in 'results.json' in the output directory. This file can later be read again
# if necessary to perform the subsequent analysis.
results = c.run(
    misfit_types=list(THRESHOLDS.keys()),
    output_directory=output_directory,
    # Setting this to True will create A LOT of debug plots which are very useful
    # to tune the algorithm to the problem at hand. For production runs set this to
    # False.
    save_debug_plots=True,
    output_format=output_format)

# This produces all kinds of plots for all components and misfits it encounters.
results.plot_all(
    output_directory=output_directory,
    thresholds=THRESHOLDS,
    event_file=event_file,
    output_format=output_format)
