function p05_reco_loop
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
para(1).scan_path = '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_41';    
para(1).bin = 4;
para(1).out_path = '/gpfs/petra3/scratch/moosmanj';

para(2).scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1400';
para(2).bin = 4;
para(2).rot_axis_offset = 19.5;
para(2).out_path = '/gpfs/petra3/scratch/moosmanj';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over all parameter sets
for nn = 1:numel( para )    
    external_parameter = para(nn);    
    p05_reco    
end