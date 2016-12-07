% Script to loop over a given list of data sets that shall be
% reconstructed.
%
% Writtten by J. Moosmann, 2016-12-07

% Cell of strings of paths to raw data sets. If string is not unique is
% interpreted as a pattern and the loop will run for all data sets which
% match the pattern.
samples = {
'/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_01'
'/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_02'
'/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_03'
'/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_05'
'/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_straw_0'
'/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_straw_2_0'
};

out_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/ImmersionTests/';

% Create cell containing all paths to loop over
samp_path = {};
for nn = 1:numel( samples )
    p =  dir( [samples{nn} '*']);
    for mm = 1:numel( p )
        samp_path = cat(2, samp_path, {sprintf('%s/%s', p(mm).folder, p(mm).name )});
    end    
end

fprintf( 'Samples to be reconstructed:\n' ) 
fprintf( ' %s\n', samp_path{:})
fprintf( 'Path to reconstructions:\n %s\n', out_path )
