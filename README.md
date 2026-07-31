# matlab
MATLAB scripts and functions

Installation of reconstruction routine for micro-tomography data
acquired at the P05 imaging beamline (IBL) or P07 high energy material
science beamline (HEMS) at PETRA III at DESY, both operated by
Helmholtz-Zentrum Hereon.

1) Log in to a GPU node on the HPC cluster MAXWELL at DESY. Nodes in
   the partitons 'max-fs-display' and 'max-display' can be accessed from
   anywhere. For more information on how to connect to the cluster, see
   https://docs.desy.de/maxwell/documentation/fastx4/. However these
   display nodes come with a load manager which kills your processes
   if they consume too much memory of CPU power.

   If this happens, you have to allocate and connect to a node on a
   different partition which is more powerful. First connect to the
   MAXWELL cluster (via a display node or similar). In a terminal
   execute  the following command to allocate a  node on our Hereon
   partition 'hzg' non-exclusively (oversubscribe): 

   salloc --partition=hzg --oversubscribe --time=7-00:00:00

   In case you don't get an allocation on our 'hzg' partition you can
   include other paritions with more nodes but which can be less
   powerful. Modify the parameter if necessary.

   salloc --partition=hzg,allgpu --mem=500G -c 40 --oversubscribe --time=7-00:00:00

   Then connect to the allocated node using ssh:

   ssh -Y $SLURM_NODELIST


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


3) How to start the tomographic reconstruction routine:

   (3.1) Change directory to the 'matlab' folder:
   
   cd matlab

   (3.2) Start MATLAB using the following script:

   ./startmatlab.sh

    The script automatically sets environment variables in order to
    use a local installation of the ASTRA toolbox, starts MATLAB, and
    sets the MATLAB search paths. (Note that this will overwrite local
    MATLAB user settings, for details see 'startup.m'.) 

    (If 'startmatlab.sh' is not executable: chmod +x startmatlab.sh)

    (3.3) If not already open, open 'p05_reco' (also for P07 scans) in
    MATLAB    located at './matlab/experiments/p05/' e.g. enter 'edit
    p05_reco' in  MATLAB's command line.   

   (3.4) Edit/check the reconstruction parameters, at least you have to
   modify 'par.scan_path'.

   (3.5) Start the reconstruction by one of the following:
    - Type 'p05_reco' in MATLAB's command line
    - Click 'RUN' button in editor tab
    - Type 'F5' key when cursor focus is in the 'p05_reco.m' file
    

    (3.6) Optinal: The reconstruction can be automatically looped over
    all data sets acquired during a beamtime and/or over different
    reconstruction parameters. How to set up a loop script to
    reconstruct several data sets:
    help p05_create_reco_loop
