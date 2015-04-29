%% pmedf_findMnePos --- return word nb of the given name in the string 
%%
%% Usage:
%%      n = pmedf_findMnePos ( motorscounters_mne_string, name )
%%
%% Given the string of 'counter_mne' or 'motor_mne' (from EDF header), it
%% returns word position number of the counter or motor named 'name'. They
%% are enumerated from 1.
%% If such a name does not exist there, then -1 is returned.
%% Technically, given a string with a 'name' as the n-th word, then this
%% routine returns n.
%% 
%% Example 1:
%%	n = pmedf_findMnePos('x y omega rot', 'omega');
%% Example 2:
%%	[h, a] = pmedf_read('scan_0000.edf');
%%	counters.mne = pmedf_findInHeader(h, 'counter_mne', 'string');
%%      monitor.pos = pmedf_findMnePos(counters.mne, 'mon');
%%
%% Author: Petr Mikulik
%% Version: 11. 5. 2008

function omega_mnepos = pmedf_findMnePos ( motors_mne, motor_name )

% Technotes:
% Don't use strtok as it is slow due to implementation in m, and also
% don't use any other loop. Thus, let us base our search on find()ing
% spaces. Suppose we seach for motor "omega". Then, we will search for
% " omega " instead.

mne2 = [' ', motors_mne, ' ']; % prepend and append a space
motor_name = [' ', motor_name, ' ']; % prepend and append a space
i = strfind(mne2, motor_name); % index of word " omega "
if i==0 % motor name not found
    omega_mnepos = -1;
    return
end

k = find(mne2==' '); % where are spaces?
omega_mnepos = find(i==k); % the given motor is the m'th word

return

% eof pmedf_findMnePos.m
