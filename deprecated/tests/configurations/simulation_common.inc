begin SETTINGS
   use_run_mode = FORWARD_MODEL
   use_lsi      = False
end SETTINGS

# Use different spectral windows than default
include = <CONTROL_FILE_DIRECTORY>/simulation_spectral_windows.inc
include = <modifiers_dir>/inp/set_spectral_windows.inc

# Use scenario as truth
include = <CONTROL_FILE_DIRECTORY>/scenario_truth_common.inc

# Remove any retrieval keywords that will mess up forward model mode when not resampling
include = <modifiers_dir>/inp/remove_retrieval_keywords.inc
