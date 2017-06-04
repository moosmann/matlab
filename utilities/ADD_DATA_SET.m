function ADD_DATA_SET( as, parameter_cell_array )
% Add (almost) all parameters, i.e. MATLAB variables given in the current
% workspace at the moment of calling this function to the cell array struct
% 'PARAMETER_CELL'.
%
% ARGUMENTS
% as : string. default: ''. restore parameter values
% back to default values after a data set is added to the loop.
% parameter_cell_array : cell array of the names of the varibles to
% used to create the parameter struct.
%
% Written by Julian Moosmann, 2017-06-02, last mod: 2017-06-03
%
% ADD_DATA_SET( as, parameter_cell_array )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    as = 0;
end
if nargin < 2
    parameter_cell_array = evalin( 'caller', 'who');
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global PARAMETER_CELL
global DEFAULT
global DATA_SET_NUM

% Create struct from variables currently available in caller workspace
for nn = 1:numel( parameter_cell_array )
    name = parameter_cell_array{nn};
    
    if ~sum( strcmpi( name, {'PARAMETER_CELL', 'DEFAULT', 'ans', 'PARAMETER_STRUCT', 'raw', 'SUBSETS', 'RUN_RECO', 'PRINT_PARAMETERS'} ) )
        parameter_struct.(name) = evalin( 'caller', name );
    end
end


if isempty( DATA_SET_NUM )
    %% Initialize counter variable and define default
    DATA_SET_NUM = 0;
    DEFAULT = parameter_struct;
else    
    %% Add parameter set to cell array of parameter structs
    DATA_SET_NUM = DATA_SET_NUM + 1;
    PARAMETER_CELL{DATA_SET_NUM} = parameter_struct;
end

%% Restore default parameters
if as == 1
    % Get default struct from caller
    default = evalin( 'caller', 'DEFAULT' );
    % Set parameter in caller back to defaults
    defn = fieldnames( default );
    for nn = 1:numel(defn)
        name = defn{nn};
        if ~strcmpi( name, 'DATA_SET_NUM' )
            assignin('caller', name, default.(name))
        end
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
