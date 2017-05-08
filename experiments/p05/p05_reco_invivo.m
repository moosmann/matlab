function p05_reco_invivo(nums)
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
% Written by Julian Moosmann. First version: 2017-02-15. Last: 2017-04-24

if nargin < 1
    nums = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sets to loop over %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw = '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/';
nn = 0;
default.interactive_determination_of_rot_axis = 0;
default.visualOutput = 0;
default.scan_path = '';
default.raw_bin = 1;
default.excentric_rot_axis = 0;
default.crop_at_rot_axis = 0;
default.stitch_projections = 0; 
default.proj_range = 1; 
default.ref_range = 1; 
default.correlation_method =  'diff';
default.do_phase_retrieval = 1;
default.phase_retrieval_method = 'tie';'qp';'qpcut';
default.phase_retrieval_reg_par = 2.5; 
default.phase_retrieval_bin_filt = 0.15; 
default.phase_retrieval_cutoff_frequ = 1 * pi; 
default.phase_padding = 1; 
default.do_tomo = 1;
default.ring_filter = 1;
default.rot_axis_offset = [];
default.rot_axis_tilt = [];
default.parfolder = '';
default.write_to_scratch = 0;
default.write_flatcor = 1;
default.write_phase_map = 0; 
default.write_sino = 0; 
default.write_sino_phase = 0; 
default.write_reco = 1; 
default.write_float = 1; 
default.write_float_binned = 1; 
default.write_8bit = 1;
default.write_8bit_binned = 1;
default.write_16bit = 0; 

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_23_00'; % rot_axis_offset 5.75;
'/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_23_01'; % Nothobranchius furzeri
'/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_23_00'; % Nothobranchius furzeri; bewegung
'/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_13_00'; % no conspicuous movement artifacts, but cell shape are unclear and nuclei not visible
'/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_13_09'; % dead
'/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_07_00'; % strong movement
'/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_13_00';

nn = nn + 1;
para(nn) = default;
para(nn).scan_path = [raw 'szeb_41'];
para(nn).rot_axis_offset = 0 / para(nn).raw_bin;
para(nn).rot_axis_tilt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nDATA SETS:')
for nn = 1:numel( para)
    [~, name] = fileparts( para(nn).scan_path );
    fprintf('\n%3u : %s', nn, name )
end

% Loop over parameter sets
if nums ~= 0
    fprintf( '\n\nSTART LOOPING \n')
    for nn = 1:numel( nums )
        num = nums(nn);
        
        external_parameter = para(num);
        [~, name] = fileparts( external_parameter.scan_path );
        fprintf('\nRECONSTRUCTION OF DATA SET NUMBER : %u : %s\n', num, name )
        p05_reco
        
    end
    fprintf( '\nRECONSTRUCTION LOOP FINISHED')
end
fprintf('\n')
