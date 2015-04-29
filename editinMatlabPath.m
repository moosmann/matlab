function editinMatlabPath(Filename,SubFolder,MatlabPath)
% Create or edit a file 'Filename' in the users MATLAB search path.

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    SubFolder = '';
end
if nargin < 3
    %MatlabPath = '/mnt/tomoraid-LSDF/users/moosmann/matlab';
    MatlabPath = [getenv('HOME') '/matlab'];
end
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Append extension
if Filename(end) ~= 'm'
    Filename = [Filename '.m'];
end
% Store present working directory
CurrentPath = pwd;
% Create and open file
eval(sprintf('cd %s/%s',MatlabPath,SubFolder));
eval(sprintf('edit %s',Filename));
eval(sprintf('cd %s',CurrentPath));
