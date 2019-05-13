function vol = stitch_volumes( scan_path, reco_subfolder )
% Stitch reconstructed volumes using log file information. Note: algorithm
% is not memory efficient. At least twice the memory of the volumes to be
% stitched is required.
%
% scan_path: string. absolute path pattern common for all height levels
% reco_subfolder: string. subfolder to 'reco' to be used for stitching.
%   E.g. 'float_rawBin2'
%
% Written by J. Moosmann
%
% vol = stitch_volumes( scan_path, reco_subfolder )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    %scan_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn22_77L_Mg5Gd_8w';% upward
    scan_path = '/asap3/petra3/gpfs/p05/2018/data/11004263/processed/syn004_96R_Mg5Gd_8w'; % downward
end
if nargin < 2
    reco_subfolder = 'float_rawBin2';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = find_volume_stitch_parameter(  scan_path, reco_subfolder, 1 );

vol = cat(3, s(s(1).stitch_order).vol );

