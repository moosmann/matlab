function synchroload2019july_11006704( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Created on 05-Aug-2019 by moosmanj

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2019/data/11006704/raw/';
visual_output = 0;
interactive_mode.rot_axis_pos = 0;

scan_path = [raw_path 'syn001_35R_Ti_8w_000']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_001']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_002']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_003']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_004']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_005']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_006']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_007']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_008']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_009']; ADD

scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_000']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_001']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_002']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_003']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_004']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_005']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_006']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_007']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_008']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_009']; ADD

tomo.rot_axis.offset = 4  / raw_bin * 0.950000;
scan_path = [raw_path 'syn003_47R_Ti_12w_000']; ADD
tomo.rot_axis.offset = 4  / raw_bin * 1.050000;
scan_path = [raw_path 'syn003_47R_Ti_12w_001']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_002']; ADD
tomo.rot_axis.offset = 4  / raw_bin * 1.450000;
scan_path = [raw_path 'syn003_47R_Ti_12w_003']; ADD
tomo.rot_axis.offset = 4  / raw_bin * 0.850000;
scan_path = [raw_path 'syn003_47R_Ti_12w_004']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_005']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_006']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_008']; ADD
tomo.rot_axis.offset = 4  / raw_bin * 1.210000;
scan_path = [raw_path 'syn003_47R_Ti_12w_012']; ADD
tomo.rot_axis.offset = 4  / raw_bin * 1.400000;
scan_path = [raw_path 'syn003_47R_Ti_12w_013']; ADD

ref_range = 6; 
raw_bin = 5;
tomo.rot_axis.offset = 1.1 * 4  / raw_bin;
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_000']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_001']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_002']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_003']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_004']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_005']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_006']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_007']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_008']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_009']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_010']; ADD
%scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_011_old']; incomplete scan
tomo.rot_axis.offset = 1.2 * 4  / raw_bin;
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_011']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_012']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_013']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_014']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_015']; ADD

tomo.rot_axis.offset = 1.0 * 4  / raw_bin;
scan_path = [raw_path 'syn005_40L_Peek_12w_000']; ADD

tomo.rot_axis.offset = 0.62 * 4  / raw_bin;
scan_path = [raw_path 'syn007_47L_Peek_12w_000']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_001']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_002']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_003']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_004']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_005']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_006']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_007']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_008']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_009']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_010']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_011']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_012']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_013']; ADD

interactive_mode.rot_axis_pos = 1;
tomo.rot_axis.offset = 0.5 * 4  / raw_bin;
scan_path = [raw_path 'syn008_47R_Ti_12w_001']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_002']; ADD
interactive_mode.rot_axis_pos = 0;
scan_path = [raw_path 'syn008_47R_Ti_12w_003']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_004']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_005']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_006']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_007']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_008']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_009']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_010']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_011']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_012']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_013']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_014']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_015']; ADD
%scan_path = [raw_path 'syn008_47R_Ti_12w_016']; ADD % Problems with scan
%log
scan_path = [raw_path 'syn008_47R_Ti_12w_017']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_018']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
