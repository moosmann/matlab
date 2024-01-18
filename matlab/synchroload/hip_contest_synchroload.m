clear all

%% Read displacement field from 3 3D volume tif stacks
dirbt = '/asap3/petra3/gpfs/p05/2021/data/11008741/processed/';
%dirsam = 'syn0117_40R_Ti_12w_000/optflow/ref_117_001_lg3_g1';
dirsam = 'syn0128_67L_Mg10Gd_8w_000/optflow/ref_128_001_lg3_g1';
p = [dirbt dirsam];
fprintf('\nsample path:\n %s',p)

tmp = regexprep(dirbt,{'/asap3','/petra3','/gpfs','/processed','/data'},'');
namesam = regexprep([tmp(2:end) dirsam],{'/'},'_');
fprintf('\nsample identifier:\n %s',namesam)

fprintf('\nReading displacement field')
v(3,:,:,:) = read_images_to_stack([p '/dz/'],1,'*.tif',[],1,1);
v(2,:,:,:) = read_images_to_stack([p '/dy/'],1,'*.tif',[],1,1);
v(1,:,:,:) = read_images_to_stack([p '/dx/'],1,'*.tif',[],1,1);

%% Export to Azizo 3D vector field
bgfs = '/beegfs/desy/user/moosmanj/hip_contest_synchroload/';
outpath = [bgfs namesam];
CheckAndMakePath(outpath)
dtype = '';
bbox = [];

fn = [outpath '/vec.am'];
fprintf('\nAvizo vector field:/n %s',fn)
fprintf('\nWriting Avizo vector field')
write_HxUniformVectorField3(fn,v,dtype,bbox,1);

%% Link warped stack to beegfs
linkname = [outpath '/warped'];
targetwarped = [p '/warped'];
s = sprintf('ln -s -T %s %s',targetwarped,linkname);
fprintf('\nLinking warped directory:\n %s\n',s)
unix(s)

%% Read convex hull for masking
pch = [dirbt '/syn0128_67L_Mg10Gd_8w_000/convex_hull/'];
m = read_images_to_stack(pch,1,'*.tif',[],1,1);

%% Warped int8
w = read_images_to_stack(targetwarped,1,'*.tif',[],1,1);
%wuint8 = FilterOutlier(w.*m,0.005,'uint8',0,1);
wuint8 = 2^8 * normat(w);

outpath2 = [outpath '/warped_uint8'];
CheckAndMakePath(p)
parfor n = 1:size(w,3)
    fn = sprintf('%s/warped_%04u.tif',outpath2,n);
    imwrite(wuint8(:,:,n),fn,'tif');
end


