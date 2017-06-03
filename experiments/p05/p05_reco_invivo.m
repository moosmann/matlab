function p05_reco_invivo( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
if nargin < 1
    SUBSETS = [];
end
if nargin < 2
    RUN_RECO = 0;
end
if nargin < 3
    PRINT_PARAMETERS = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interactive_determination_of_rot_axis = 0;
visualOutput = 0;
scan_path = '';
raw_bin = 1;
excentric_rot_axis = 0;
crop_at_rot_axis = 0;
stitch_projections = 0; 
proj_range = 1; 
ref_range = 1; 
correlation_method =  'diff';
do_phase_retrieval = 1;
phase_retrieval_method = 'tie';'qp';'qpcut';
phase_retrieval_reg_par = 2.5; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
phase_padding = 1; 
do_tomo = 1;
ring_filter = 1;
rot_axis_offset = [];
rot_axis_tilt = [];
parfolder = '';
write_to_scratch = 0;
write_flatcor = 1;
write_phase_map = 0; 
write_sino = 0; 
write_sino_phase = 0; 
write_reco = 1; 
write_float = 1; 
write_float_binned = 1; 
write_8bit = 1;
write_8bit_binned = 1;
write_16bit = 0; 
ADD_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_path = '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/';

'szeb_23_00'; % rot_axis_offset 5.75;
'szeb_23_01'; % Nothobranchius furzeri
'szeb_23_00'; % Nothobranchius furzeri; bewegung
'szeb_13_00'; % no conspicuous movement artifacts, but cell shape are unclear and nuclei not visible
'szeb_13_09'; % dead
'szeb_07_00'; % strong movement
'szeb_13_00';


scan_path = [raw_path 'szeb_41'];
rot_axis_offset = 0 / raw_bin;
rot_axis_tilt = 0;
ADD_DATA_SET

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS, PARAMETER_CELL)
