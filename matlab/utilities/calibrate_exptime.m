scan_path = '/asap3/petra3/gpfs/p07/2021/data/11009678/raw/testfly012';
%scan_path = pwd;
[par_path, name] = fileparts( scan_path);
nx = sprintf( '%s/%s_nexus.h5', scan_path, name );

n_angle = h5read( nx,'/entry/scan/n_angle');
n_dark = h5read( nx,'/entry/scan/n_dark');
n_ref = h5read( nx,'/entry/scan/n_ref');
n_img = h5read( nx,'/entry/scan/n_img');
exptime = h5read( nx,'/entry/hardware/camera1/exptime') / 1000;
rotation = h5read( nx,'/entry/scan/mode');
srt = double( h5read( nx, '/entry/scan/data/s_rot/time' ) );
srv = double( h5read( nx, '/entry/scan/data/s_rot/value' ) );
imk = h5read( nx,'/entry/scan/data/image_key/value');
imk_srot = imk(1+n_dark:end);

srot_time = srt(imk_srot == 0) / 1000;
srot_val = srv(imk_srot == 0);

total_time = srot_time(end) - srot_time(1);
total_angle = srot_val(end) - srot_val(1);
scalfac = 1.0001;
total_time_set = n_angle * exptime * scalfac;
total_angle_set = rotation;

o = total_time / n_angle / scalfac - exptime;

fprintf( '\nn_angle: %u', n_angle )
fprintf( '\nn_dark: %u', n_dark)
fprintf( '\nn_ref: %u', n_ref)
fprintf( '\nexptime: %f', exptime * 1000)
fprintf( '\nnumel(s_rot_time): %u', numel(srt))
fprintf( '\nnumel(imk): %u', numel(imk))
fprintf( '\nnumel(imk_srot): %u', numel(imk_srot))

fprintf( '\ntotal time: %f (is), %f (set)', total_time, total_time_set )
fprintf( '\ntotal angle: %f (is), %f (set), %f (set no overlap)(', total_angle, total_angle_set, total_angle_set -  total_angle_set / n_angle)

fprintf( '\noffset: %f ms', o * 1000)




fprintf( '\n' )