function p05_reco_invivo(nums, doreco)
if nargin < 1
    nums = [];
end
if nargin < 2
    doreco = 0;
end
if nargin < 3
    print_field = '';
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
p05_reco_loop( nums, doreco, print_field, para)
