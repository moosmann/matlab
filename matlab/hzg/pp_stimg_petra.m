function [stimg_name, stimg_key, petra, petra_scan] = pp_stimg_petra(nexuslog_name,par)
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
    stimg_name_value_n = unique( h5read( nexuslog_name{n}, '/entry/scan/data/image_file/value') );
    stimg_name_value = cat( 1, stimg_name_value, stimg_name_value_n);

    stimg_name_time_n = h5read( nexuslog_name{n},'/entry/scan/data/image_file/time');
    stimg_name_time = cat( 1,stimg_name_time, stimg_name_time_n);

    stimg_key_value_n = h5read( nexuslog_name{n},'/entry/scan/data/image_key/value');
    stimg_key_value = cat( 1,stimg_key_value,stimg_key_value_n);

    stimg_key_time_n = h5read( nexuslog_name{n},'/entry/scan/data/image_key/time');
    stimg_key_time = cat( 1,stimg_key_time, stimg_key_time_n);

    % PETRA ring current
    if par.ring_current_normalization == 1
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
    end
    %if isscalar(nexuslog_name)
    stimg_name.scan(n).value = stimg_name_value_n;
    stimg_name.scan(n).time = stimg_name_time_n;
    stimg_key.scan(n).value = stimg_key_value_n;
    stimg_key.scan(n).time = stimg_key_time_n;

    if par.ring_current_normalization == 1
        petra_scan.time = pt;
        petra_scan.current = pc;
    end
    %end
end

stimg_name.time = stimg_name_time;
stimg_name.value = stimg_name_value;
stimg_key.time = stimg_key_time;
stimg_key.value = stimg_key_value;
% Sort PETRA time
if par.ring_current_normalization == 1
    [petra.time, index]  = sort( petra_time );
    petra.current = petra_current(index);
else
    petra = [];
end
