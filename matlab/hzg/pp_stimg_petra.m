function [stimg_name, stimg_key, petra, petra_scan] = pp_stimg_petra( nexuslog_name )
% Helper function returning name, time stamp and image keys as well as
% petra III current from the nexus log.

% Initialize
stimg_name_time = [];
stimg_name_value = {};
stimg_key_time = [];
stimg_key_value = [];
petra_time = [];
petra_current = [];
petra_scan = [];

% Loop over nexus log files
for n = 1:numel( nexuslog_name )
    % Images name and key
    stimg_name_value = cat( 1, stimg_name_value, unique( h5read( nexuslog_name{n}, '/entry/scan/data/image_file/value') ) );
    stimg_name_time = cat( 1, stimg_name_time, h5read( nexuslog_name{n},'/entry/scan/data/image_file/time') );
    stimg_key_value = cat( 1, stimg_key_value, h5read( nexuslog_name{n},'/entry/scan/data/image_key/value') );
    stimg_key_time = cat( 1, stimg_key_time, h5read( nexuslog_name{n},'/entry/scan/data/image_key/time') );
    
    % PETRA ring current
    [pt, index] = unique( h5read( nexuslog_name{n},'/entry/hardware/beam_current/current/time') );
    pc = h5read( nexuslog_name{n},'/entry/hardware/beam_current/current/value');
    pc = pc(index);
    nn = 1;
    while pc(nn) == 0 || pt(nn) == 0
        pt(nn) = [];
        pc(nn) = [];
        nn = nn + 1;
    end        
    petra_time = cat( 1, petra_time, pt );
    petra_current = cat( 1, petra_current, pc );
    if n == 1
        stimg_name.scan.time = stimg_name_time;
        stimg_name.scan.value = stimg_name_value;
        stimg_key.scan.time = stimg_key_time;
        stimg_key.scan.value = stimg_key_value;
        petra_scan.time = pt;
        petra_scan.current = pc;
    end
end

stimg_name.time = stimg_name_time;
stimg_name.value = stimg_name_value;
stimg_key.time = stimg_key_time;
stimg_key.value = stimg_key_value;
% Sort PETRA time
[petra.time, index]  = sort( petra_time );
petra.current = petra_current(index);
