function editinMatlabPath(filename,subfolder,MatlabPath)
% Create or edit a file 'filename' in the users MATLAB search path.
%
% Written by Julian Moosmann. Last modification: 2016-09-28
%
% editinMatlabPath(filename,subfolder,MatlabPath)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    subfolder = '';
end
if nargin < 3
    %MatlabPath = '/mnt/tomoraid-LSDF/users/moosmann/matlab';
    MatlabPath = [getenv('HOME'), filesep, 'matlab'];
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Append extension
if filename(end) ~= 'm'
    filename = [filename '.m'];
end
% Store present working directory
CurrentPath = pwd;
% Create and open file
eval( sprintf('cd %s%s%s', MatlabPath, filesep, subfolder) );
eval( sprintf('edit %s', filename) );
eval( sprintf('cd %s', CurrentPath) );
