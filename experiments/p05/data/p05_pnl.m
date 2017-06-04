%% Script to check P05 data on sponges. Problem: FOV extension was used, but rotation axis excentricity was beyond FOV.

scan = 'pnl_01_e';

% parent directory (raw)
par_dir = '/asap3/petra3/gpfs/p05/2016/data/11001464/raw/';

% data directory
data_dir = [par_dir scan '/'];

cd(data_dir)

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

% Dark field
verbose = 0;
read_dark = @(num) FilterPixel( read_dat( sprintf('%s%s', data_dir, dark_names{num}) ), [0.02 0.02], verbose);
for nn = num_dark:-1:1
    stack(:, :, nn) = read_dark(nn);
end
dark = FilterPixel( squeeze( median(stack, 3) )', [0.02 0.02], verbose);
clear stack;

% Flat field
read_ref = @(num) FilterPixel( read_dat( sprintf('%s%s', data_dir, ref_names{num}) ) )';
for nn = num_ref:-1:1
    stack(:, :, nn) = read_ref(nn);
end


read_img = @(num) FilterPixel( read_dat( sprintf('%s%s', data_dir, img_names{num}) ) )';



p1 = read_img(1);
p2 = read_img( floor(num_img / 2) );
f1 = read_ref(1);
f2 = read_ref( floor(num_ref / 2) );

% print info
domain(dark)
domain(f1)
domain(f2)
domain(p1)
domain(p2)

% flat-dark correction, roi
x = 231:880;
make_im = @(proj,dark,flat) (proj(x,:) - dark(x,:)) ./ (flat(x,:) - dark(x,:));

im1 = make_im(p1, dark, f1);
im2 = make_im(p2, dark, f2);

itool(im1)
itool(fliplr(im2))