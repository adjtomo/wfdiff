# To run this script
# mpirun -n 8 python run_wfdiff_test.py

from wfdiff.wfdiff import WFDiff

# A non-existant directory where the output will be stored.
#OUTPUT_DIRECTORY = "output"

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


# XXX : Move this specfem_util.py
# SPECFEM encoded information in the filenames and the way this
# is done differs depending on the SPECFEM version. This function
# will be used to extract network id, station id, and the component
# from any waveform file encountered in case SPECFEM ASCII output
# is used.
def get_net_sta_comp(filename):
    net, sta, chan, _ = filename.split(".")
    return net, sta, chan[-1]

#--------------------------------------------------------------------------------
# Run example

iex = 1

if iex == 1:
    # test data (using ENZ components)
    low_res_seismos="./test_data/NGLL5/*semd"
    high_res_seismos="./test_data/NGLL7/*semd"
    station_info="./test_data/STATIONS"
    cmtsolution_file='./test_data/CMTSOLUTION'
    trace_tags = ['NGLL5','NGLL7']
    is_specfem_ascii = True
    OUTPUT_DIRECTORY = "output_test"
    outformat = 'pdf' # or ps, eps etc.

elif iex == 2: 
    # Test using sac data (using ENZ components) - Same as the Example 1
    # Run specfem_util_test.py first to create rotate SAC file for NGLL5 and NGLL7 case
    # test (using ENZ components) - Same as the Example 1
    low_res_seismos="./20140831030657110_NGLL5_ENZ/*sac"
    high_res_seismos="./20140831030657110_NGLL7_ENZ/*sac"
    trace_tags = ['NGLL5-ENZ','NGLL7-ENZ']
    is_specfem_ascii = False
    OUTPUT_DIRECTORY = "output_test"
    outformat = 'pdf' # or ps, eps etc.

elif iex == 3:
    # Test using sac data (using RTZ components)
    # Run specfem_util_test.py first to create rotate SAC file for NGLL5 and NGLL7 case
    low_res_seismos="./20140831030657110_NGLL5_RTZ/*sac"
    high_res_seismos="./20140831030657110_NGLL7_RTZ/*sac"
    station_info="./test_data/STATIONS"
    trace_tags = ['NGLL5-RTZ','NGLL7-RTZ']
    is_specfem_ascii = False
    OUTPUT_DIRECTORY = "output_test"
    outformat = 'pdf' # or ps, eps etc.

#--------------------------------------------------------------------------------

# Configuration.
c = WFDiff(
    low_res_seismos=low_res_seismos, high_res_seismos=high_res_seismos, station_info=station_info,
 
    # Specify the units of the data and the units the analysis should take
    # place in.
    data_units="displacement", desired_analysis_units="displacement",
    # Periods to test.
    t_min=1, t_max=10, dt=1,
    get_net_sta_comp_fct=get_net_sta_comp,
    # Set to True if waveform ASCII files are used. All other fileformat
    # should otherwise work just fine.
    is_specfem_ascii=is_specfem_ascii,
    # Data window to take into account in seconds since the first sample.
    starttime=0, endtime=100)

# Perform the frequency dependent misfit measurements. Results will be stored
# in 'results.json' in the output directory. This file can later be read again
# if necessary to perform the subsequent analysis.
results = c.run(
    misfit_types=list(THRESHOLDS.keys()),
    output_directory=OUTPUT_DIRECTORY,trace_tags=trace_tags,
    # Setting this to True will create A LOT of debug plots which are very useful
    # to tune the algorithm to the problem at hand. For production runs set this to
    # False.
    save_debug_plots=True,
    outformat=outformat)

# This produces all kinds of plots for all components and misfits it encounters.
results.plot_all(
    output_directory=OUTPUT_DIRECTORY,
    thresholds=THRESHOLDS,
    cmtsolution_file=cmtsolution_file,
    outformat=outformat)
