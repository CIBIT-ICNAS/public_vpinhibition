# scripts-02-stimulus
MATLAB code for generating the stimuli for the experiment.

Dependencies:
- Psychtoolbox 3.0.14
- Neuroelf v1.1 rc2

## Files:
- `MAIN_Behavioral.m` - Main script for generating the stimuli for the behavioral experiment.
- `MAIN_MRI.m` - Main script for generating the stimuli for the fMRI experiment.

## Folders:
- `functions` - helper functions used by the main functions.
- `functions_main` - main stimulus implementation.
- `functions_mri` - specific functions for the MRI setup (eyetracker and trigger).
- `input` - input files for the experiment (protocol for each run and texture with the plaid)
- `output` - output files for the experiment (log files and saved workspace)
- `prt` - protocol (PRT) files for each run in BrainVoyager format (output of `scripts-01-protocolcreator`).

## How to:
1. Open Main_Behavioral.m
    1. Change subject name and initials (template <S00>)
    2. Run section by section
    3. In the end, run the "Save Workspace" section, for joint output of the keypresses.

## Tested on:
* Windows 10 Pro Version 1809
* MATLAB R2017b
* NeuroElf v1.1 rc2
* Psychtoolbox 3.0.14
