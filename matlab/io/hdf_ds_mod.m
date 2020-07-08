clear all

fn_org =  '/asap3/petra3/gpfs/p05/2020/data/11008476/raw/hzb_021_GH38powder20/hzb_021_GH38powder20_nexus.h5';
fn_fix =  '/asap3/petra3/gpfs/p05/2020/data/11008476/raw/hzb_021_GH38powder20/hzb_021_GH38powder20_nexus_fix.h5';

%exist( fn_org, 'file' )
%exist( fn_fix, 'file' )

%% Read old table
val = h5read( fn_org, '/entry/scan/data/image_file/value' );
h5disp( fn_org, '/entry/scan/data/image_file/value' );

%% Extent table
val2(2:numel(val) + 1,1) = val;
val2(1,1) = "hzb_021_GH38powder20_000000_dar.tif";

%% Extent data set in HDF%

fid = H5F.open( fn_fix,'H5F_ACC_RDWR','H5P_DEFAULT');
dset_id = H5D.open(fid, '/entry/scan/data/image_file/value' );
H5D.set_extent(dset_id,2551);
H5D.close(dset_id);
H5F.close(fid);

%% Write extented data set to hdf
h5write( fn_fix, '/entry/scan/data/image_file/value', val2 );