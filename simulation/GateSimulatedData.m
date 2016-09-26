%% Data sets

% sino matrix: the object is initialized with the contents of this matrix.
% The matrix must be of size (u,angles,v), where u is the number of columns
% of the detector and v the number of rows as defined in the projection
% geometry 

data(1).dir = '/home/jmoosmann/data/gate/';
data(1).filename = 'detector_two_spheres_Astra_20150313';
data(1).fieldname = 'detector_astra';
data(1).permuteOrder = [3 2 1];
data(1).fullAngle_rad = pi;

data(2).dir = '/home/jmoosmann/data/gate/';
data(2).filename = '20150317_water_spheres_all_photons';
data(2).fieldname = 'detector_astra';
data(2).permuteOrder = [1 2 3];
data(2).fullAngle_rad = 2 * pi;

data(3).dir = '/home/jmoosmann/data/gate/';
data(3).filename = '20150317_water_spheres_primary_photons';
data(3).fieldname = 'detector_astra';
data(3).permuteOrder = [1 2 3];
data(3).fullAngle_rad = 2 * pi;

data(4).dir = '/home/jmoosmann/data/gate/';
data(4).filename = '20150317_water_spheres_all_photons_high_act';
data(4).fieldname = 'detector_astra';
data(4).permuteOrder = [1 2 3];
data(4).fullAngle_rad = 2 * pi;

