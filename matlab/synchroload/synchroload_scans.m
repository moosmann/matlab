%% Beamtimes
beamtime = { ...
    2016 11001978
    2017 11003950
    2017 11004016
    2017 11003288
    2017 11003440
    2017 11003773
    2018 11004263
    2018 11004936
    2018 11005553
    2019 11005842
    2019 11006704
    2019 11006991
    2020 11010056
    2020 11010206
    2021 11009431
    2021 11008741
    };

%% Output directory
%segpath = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/segmentation/Ti';

%% Read scans %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over beamtimes
s(1).name = '';
cc = 1;
for nn = 1:size( beamtime, 1 )
    
    % beamtime folder
    if size( beamtime, 1 ) == 1
        beamtime_year = beamtime{1};
        beamtime_id = beamtime{2};
    else
        beamtime_year = beamtime{nn,1};
        beamtime_id = beamtime{nn,2};
    end
    p = sprintf( '/asap3/petra3/gpfs/p05/%u/data/%u/raw', beamtime_year, beamtime_id );
    
    % scan folder
    dir_struct = dir( p );
    m = [dir_struct.isdir];
    dir_struct = dir_struct(m);
    name_pattern = {'mah','syn'};
    match = contains( {dir_struct.name}, name_pattern, 'IgnoreCase', 1 );
    scan_struct = dir_struct(match);
    
    
    % Loop over scans
    for mm = 1:length( scan_struct )
        name = scan_struct(mm).name;
        s(cc).name = name;
        s(cc).beamtime_year = beamtime_year;
        s(cc).beamtime_id = beamtime_id;
        cc = cc + 1;
        fprintf( '\n %u %u %s', beamtime_year, beamtime_id, name )
        
    end
end
