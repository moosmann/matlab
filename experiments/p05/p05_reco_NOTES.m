% Notes for P05 reconstruction pipeline 'p05_reco'.

%% Padding and phase retrieval
% Symmetric padding of intensity maps before phase retrieval cleary reduces
% artifacts in the retrieved phase maps which are due to inconistent
% boundary conditions. However, due to local tomography phase retrieval
% with padding can result in more artifacts in the recostructed volume
% (such as additonal bumps at the halo described below). After phase
% retrieval without padding the region close to the left
% boundary is blended in to the region close to the right boundary and vice
% versa (same holds for top and bottom). This effect fades out with
% increasing distance to the boundarys. Due to this blending effect
% asymmetrical 'density' distributions, which are the cause of halo-like
% artifacts stretching outwards from the biggest possible circle within the
% reconstruction volume within a slice, are reduced which reduces the
% halo-like artifacts as well.

%% FOV extension by excentric rotation axis
% For absorpion-contrast data or more precisely when no phase retrieval
% is desired, volumes can be reconstructed from a data set where an excentric
% position of the rotation axis was used without prior stitching of the
% projections by simply providing the correct rotation axis position and
% setting fbp_filter_padding to 1. Maybe, for the automatic detection of the
% rotation axis to work, the area to correlate has to be adjusted i.e.
% 'rot_corr_area1'. When using phase retrieval, this
% approach does not work appropriately and gives rise to artifacts near the
% center of the reconstructed volume. This is due to the fact, that without
% stitching the phase is retrieved from a 'cropped' projection which
% results in inconsitently retrieved low frequencies (large scale
% variations) in the phase map. Using the 'linear' FBP filter instead of
% 'Ram-Lak' maybe reduces these artifacts (not tested).

%% Correlation of projections and flat fields
% Compared to the simple differencing method, the cross correlation method
% is very likely to give non-optimal results since is more sensitive to
% small-scale features, such as those stemming from contimation of the beam
% from the scintllator, diamond window, etc, and less sensitive to
% variations on a larger scale.
% Matlab's SSIM becoms time consuming for large data sets because it
% involves Gaussian blurring.

%% Entropy-type determination of rotation axis position
% Empirically, this does not work well for excentric rotation axis positions and
% stitchted projections. To be tested

%% Ring artefact filter

%% Errors
% GPU CUDA error: device busy
% Do not call 'gpuDevice()' when ASTRA toolbox is used with multi-GPU
% support, or MATLAB terminates abnormally with a segmentation violation.
