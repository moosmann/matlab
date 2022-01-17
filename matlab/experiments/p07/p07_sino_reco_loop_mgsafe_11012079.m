
parpath = '/asap3/petra3/gpfs/p07/2021/data/11012079/processed/';

%Trans02 zur Reko bereit

%360° Rekos:
rot_angle_full_range = 2 * pi;
rot_axis_search_range = 3:0.1:5.5;
%scan_path = [parpath 'mgsafe028_544'];
%sino_reco
scan_path = [parpath 'mgsafe027_462'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
% (hattest Du wohl schon rekonstruiert)
scan_path = [parpath 'mgsafe026_542'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
% (eventuel auch schon)
scan_path = [parpath 'mgsafe025_473'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe024_589'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe023_273'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe022_272'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);

% 180° Rekos
rot_angle_full_range = 1 * pi;
rot_axis_search_range = 1:0.1:4.5;
scan_path = [parpath 'mgsafe019_541'];
%sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);

scan_path = [parpath 'mgsafe018_540'];
%sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);

scan_path = [parpath 'mgsafe015_280'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe014_588'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe012_507'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe010_503'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe009_505'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe008_447'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe006_449'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe005_279'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe004_282'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);
scan_path = [parpath 'mgsafe003_281'];
sino_reco(scan_path, rot_angle_full_range, rot_axis_search_range);