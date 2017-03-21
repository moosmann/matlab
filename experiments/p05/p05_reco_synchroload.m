function p05_reco_synchroload
% Script to loop over parameter sets related to paramters of scritp
% 'p05_reco'. Set parameters to loop over as elements of the structure
% array 'para' below. Fieldnames of 'para' MUST match the names of
% parameters in 'p05_reco'.
%
% Visual output ('visualOutput') and user interaction
% ('interactive_determination_of_rot_axis') are turned off by default if
% not set otherwise.
%
% Start loop by pushing 'F5', clicking on 'Run' in the Editor tab, or
% typing 'p05_reco_loop' in the Command Window.
%
% Written by Julian Moosmann. First version: 2017-02-15. Last: 2017-02-20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sets to loop over %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = 1;
para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_15_57R';
para(nn).bin = 2;
para(n).do_stitching = 0; 
para(nn).proj_range = 1; 
para(nn).ref_range = 1; 
para(nn).do_phase_retrieval = 0;
para(nn).do_tomo = 1;
para(nn).rot_corr_area1 = [0.25 75];

% nn = nn + 1;
% para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1400';
% para(nn).bin = 4;
% para(nn).rot_axis_offset = 19.5;
% para(nn).out_path = '/gpfs/petra3/scratch/moosmanj';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over all parameter sets
for nn = 1:numel( para )    
    external_parameter = para(nn);    
    p05_reco    
end