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

# test data (using ENZ components)
low_res_seismos="./test_data/NGLL5/*semd"
high_res_seismos="./test_data/NGLL7/*semd"
stations_file="./test_data/STATIONS"
event_file='./test_data/CMTSOLUTION'
trace_tags = ['NGLL5','NGLL7']
wf_format = 'specfem'
output_directory = "output_test"
output_format = 'pdf' # or ps, eps etc.


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
    trace_tags=trace_tags,
    # asdf dataset tags
    # asdf_tags=asdf_tags,
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
