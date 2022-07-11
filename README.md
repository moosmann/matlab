# matlab
MATLAB scripts and functions

Installation of reconstruction routine for micro-tomography data
acquired at the P05 imaging beamline (IBL) or P07 high energy material
science beamline (HEMS) at PETRA III at DESY, both operated by
Helmholtz-Zentrum Hereon.

1) Log in to a GPU node on the MAXWELL cluster at DESY,
   e.g. max-hzgg001 to max-hzgg006, max-fsg, etc. The nodes max-nova,
   max-display and max-display3 are meant for remote access and
   visualization, but not for computationally expensive
   tasks. However, very small scans or single slice can be
   reconstructed on these nodes.

2) Download the latest MATLAB files from GitHub:
   
   git clone https://github.com/moosmann/matlab.git

   OR

   Update to latest version: Fetch latest version of the 'master' branch 
   from the remote repository 'origin' and reset/overwrite all local 
   changes. New files that were created locally and which do
   not exist in the latest branch are not deleted.

   git fetch origin master & git reset --hard origin/master

   Before updating you can back-up local files by branching if necessary:

   git add --all
   git commit -m "COMMIT MESSAGE"
   git branch NAME_OF_BRANCH

   To check which files were modified and which new (untracked files):
    
   git status


How to start the tomographic reconstruction routine:

1) Change directory to the 'matlab' folder:
   
   cd matlab

2) Start MATLAB using the following script:

   ./startmatlab.sh

    The script automatically sets environment variables in order to
    use a local installation of the ASTRA toolbox, starts MATLAB, and
    sets the MATLAB search paths. (Note that this will overwrite local
    MATLAB user settings, for details see 'startup.m'.)

    (If 'startmatlab.sh' is not executable: chmod +x startmatlab.sh)

3) If not already open, open 'p05_reco' (also for P07 scans) in MATLAB
   located at './matlab/experiments/p05/' e.g. enter 'edit p05_reco' in
   MATLAB's command line.

4) Edit/check the reconstruction parameters, at least you have to
   modify 'par.scan_path'.

5) Start the reconstruction by one of the following:
    - Type 'p05_reco' in MATLAB's command line
    - Click 'RUN' button in editor tab
    - Type 'F5' key when focus is in the 'p05_reco.m' file    -


The reconstruction can be automatically looped over all data sets
acquired during a beamtime and/or over different reconstruction
parameters. How to set up a loop script to reconstruct several data
sets: help p05_create_reco_loop
