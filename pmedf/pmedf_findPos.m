%% pmedf_findPos --- return value of the counter/motor position
%%
%% Usage 1:
%%      value = pmedf_findPos ( motorscounters_pos_string, n )
%% Usage 2:
%%      value = pmedf_findPos ( header, 'motor|counter', name )
%%
%% Ad Usage 1:
%% Given the string of 'counter_pos' or 'motor_pos' (from EDF header), it
%% returns the value of the n-th counter or motor. They are enumerated from 1.
%% If the motor position is out of range, then 1e31 is returned.
%% Technically, this routine returns n-th word (converted to float) from the
%% given string.
%% 
%% Ad Usage 2:
%%  Like Usage 1, but it proceeds like Example 1 -- from the header, it gets
%% '{motor|counter}_{pos|mne}', which it parses for the given object name.
%%
%% Example 1a:
%%	omega.value = pmedf_findPos('1.23 4.56', 2);
%% Example 1b:
%%	[h, a] = pmedf_read('scan_0000.edf');
%%	% Find value of the "omega" motor:
%%	motors.mne = pmedf_findInHeader(h, 'motor_mne', 'string');
%%	motors.pos = pmedf_findInHeader(h, 'motor_pos', 'string');
%%      omega.pos = pmedf_findMnePos(motors.mne, 'omega');
%%      omega.value = pmedf_findPos(motors.pos, omega.pos);
%%	% Find value of the "mon" counter:
%%	counters.mne = pmedf_findInHeader(h, 'counter_mne', 'string');
%%	counters.pos = pmedf_findInHeader(h, 'counter_pos', 'string');
%%      monitor.pos = pmedf_findMnePos(counters.mne, 'mon');
%%      monitor.value = pmedf_findPos(counters.pos, monitor.pos);
%%
%% Example 2:
%%      omega.value = pmedf_findPos(h, 'motor', 'omega');
%%      monitor.value = pmedf_findPos(h, 'counter', 'mon');
%%
%% Author: Petr Mikulik
%% Version: 13. 2. 2005


function omega_pos = pmedf_findPos ( motors_pos, motor_mnepos, opt )

if nargin<2 | nargin>3
    error('usage: pmedf_emptyHeader( [headersize] )');
end

if nargin==3
    h = motors_pos; % EDF header
    which = motor_mnepos; % 'counter' or 'motor' prefix
    which_mne = pmedf_findInHeader(h, [which, '_mne'], 'string');
    motors_pos = pmedf_findInHeader(h, [which, '_pos'], 'string');
    motor_mnepos = pmedf_findMnePos(which_mne, opt);
end

if motor_mnepos <= 0
    omega_pos = 1e31;
    return
end

% Search for the motor_npos's word in the motor_pos string.
k = find(motors_pos==' '); % where are spaces?
% motors_pos(1:k(motor_mnepos-1)-1) = '';
if motor_mnepos>length(k)
    motors_pos = motors_pos(k(motor_mnepos-1):end);
else
    if motor_mnepos<=1
	motors_pos = motors_pos(1:k(1)-1);
    else
	motors_pos = motors_pos(k(motor_mnepos-1):k(motor_mnepos)-1);
    end
end

omega_pos = sscanf(motors_pos, '%g');

% eof pmedf_findPos.m
