%% P05 SynchroLoad data on Mg screws
% Data preprocessing so far.
% Written by Julian Moosmann.
% First version: 2016-09-28. Last modifcation: 2016-09-08.

par_dir = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160913_000_synload/raw/';
scan = 'mg10gd_41_2w';

% output folder
out_dir = [par_dir(1:end-4) 'processed/' scan '/test_jm/'];
flatcor_dir = [out_dir 'flat_corrected/'];
reco_dir = [out_dir 'reco/'];

% data directory
data_dir = [par_dir scan '/'];

cd(out_dir)

% get file names
data_struct = dir( [data_dir '*.img'] );
img_names = {data_struct.name};
num_img = numel(img_names);

data_struct = dir( [data_dir '*.ref'] );
ref_names = {data_struct.name};
num_ref = numel(ref_names);

data_struct = dir( [data_dir '*.dar'] );
dark_names = {data_struct.name};
num_dark = numel(dark_names);

%% Dark field
tic
verbose = 0;
read_dark = @(num) FilterPixel( read_dat( sprintf('%s%s', data_dir, dark_names{num}) ), [0.02 0.02], verbose);
for nn = num_dark:-1:1
    stack(:, :, nn) = read_dark(nn);
end
dark = FilterPixel( squeeze( median(stack, 3) )', [0.02 0.02], verbose);
clear stack;
fprintf('Median dark field: %g s', toc)

%% Flat field
tic
read_ref = @(num) FilterPixel( read_dat( sprintf('%s%s', data_dir, ref_names{num}) ), [0.01 0.005]);
for nn = num_ref:-1:1
    stack(:, :, nn) = read_ref(nn);
end
ref = FilterPixel( squeeze( mean(stack, 3) )', [0.005 0.0025]);
clear stack
fprintf('Mean flat field: %g s', toc)

%% Projections
tic
CheckAndMakePath(flatcor_dir)
read_img = @(num) FilterPixel( read_dat( sprintf('%s%s', data_dir, img_names{num}) ), [0.02 0.01])';
for nn = 1:num_img
    proj = ( read_img(nn) - dark ) ./ ( ref - dark );
    filename = sprintf('%sproj_%06u.tif', flatcor_dir, nn );
    write32bitTIFfromSingle(filename, proj);
end
fprintf('Save flat field corrected projections: %g s = %g min', toc, toc/60)


%% Reco
tic
CheckAndMakePath(reco_dir)
theta = 180 * (0:num_img - 1) / (num_img - 1);
s = squeeze(sino(:, :, 1))';
rec = iradon(s(6:end,:), theta);
nn = 1;
filename = sprintf('%sreco_%06u.tif', reco_dir, nn );
write32bitTIF(filename, rec)
fprintf('Save reconstructed slices: %g s = %g min', toc, toc/60)


%% Show
figure
x = 601:1400;
imshow(rec(x,x),[],'InitialMagnification','fit')

fprintf('\nFINISHED\n')