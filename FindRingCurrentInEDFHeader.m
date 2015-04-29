function srcur = FindRingCurrentInEDFHeader(Header)
% Find value of ring current if it is save in the header of an edf.
%
% !!
% Modification of function "pmedf_findPos" is needed in line 74:
% "omega_pos = sscanf(motors_pos, '%g','C')" -> "omega_pos = sscanf(motors_pos, '%g')"
% !!
%
% Written by Julian Moosmann, 2013-10-09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counters.mne = pmedf_findInHeader(Header, 'counter_mne', 'string');
if ~isempty(counters.mne)
    counters.pos = pmedf_findInHeader(Header, 'counter_pos', 'string');
    srcur_pos    = pmedf_findMnePos(counters.mne, 'srcur');
    srcur = pmedf_findPos(counters.pos, srcur_pos);
else
    srcur = [];
end