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
nn = 0;
default.scan_path = '';
default.bin = 2;
default.excentric_rot_axis = 0;
default.stitch_projections = 0; 
default.proj_range = 1; 
default.ref_range = 1; 
default.do_phase_retrieval = 0;
default.do_tomo = 1;
default.crop_at_rot_axis = 0;
default.rot_axis_offset = [];
default.rot_axis_tilt = [];
%default.rot_corr_area1 = [0.25 75];

nn = nn + 1;
para(nn) = default;
para(nn).excentric_rot_axis = 0;
para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_15_57R';

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_10_13R_top';
para(nn).excentric_rot_axis = 1;
para(nn).rot_axis_offset = 538.5;
para(nn).rot_axis_tilt = -0.00264; % about -.15 degrees

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_10_13R_bottom';
para(nn).excentric_rot_axis = 1;
para(nn).rot_axis_offset = 539;
para(nn).rot_axis_tilt = -0.00264; % about -.15 degrees

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_11_20R_top';
para(nn).excentric_rot_axis = 1;
para(nn).rot_axis_offset = 538.5;

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_11_20R_bottom';
para(nn).excentric_rot_axis = 1;
para(nn).rot_axis_offset = 538.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over all parameter sets
for nn = 1:numel( para )    
    external_parameter = para(nn);    
    p05_reco    
end