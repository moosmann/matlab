function OutputPath = MakePath(varargin)
% Define a path using 'sprintf'. It is checked if folder
% already exist and if not it is created. A trailing slash is added.
% 
% Arguments are directly passed to sprintf: formatSpec,A1,...,A2
%
% Written by Julian Moosmann, last version 2013-10-31
%
% OutputPath = MakePath(varargin)

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creat path string
OutputPath = sprintf(varargin{:});
% Check if path exists
if ~exist(OutputPath,'dir')
    mkdir(OutputPath);
end
% Add trailing slash to output variable 'OutputPath'
if OutputPath(end) ~= '/'
    OutputPath = [OutputPath '/'];
end
% Add trailing slash to input variable
if nargin == 1 && nargout == 0
    assignin('caller',inputname(1),OutputPath)
end
