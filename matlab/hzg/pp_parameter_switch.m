% List all parameters defined so far
par_quick_switch.who = who();
% Number of parameters defined
par_quick_switch.num_structs = numel( par_quick_switch.who );
% Loop of parameters
for n = 1:par_quick_switch.num_structs
    % Name of parameter struct
    sn = par_quick_switch.who{n};
    % Store parameter struct in fields
    par_quick_switch.structs.(sn) = eval( sn );
end