% Test datastore
sino_path = '/asap3/petra3/gpfs/p05/2019/data/11007580/processed/smf_09_be_3033/trans04';

ds = datastore( sino_path, ...
    'IncludeSubfolders', true, ...
    'FileExtensions', '.tif', ...
    'Type', 'image');