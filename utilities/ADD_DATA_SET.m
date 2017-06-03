function ADD_DATA_SET( restore_default_after_set, parameter_cell_array )
% Add (almost) all parameters, i.e. MATLAB variables given in the current
% workspace at the moment of calling this function to the cell struct
% 'par'.
%
% ARGUMENTS
% restore_default_after_set : bool. default: 0. restore parameter values
% back to default values after a data set is added to the loop.
% parameter_cell_array : cell array of the names of the varibles to
% used to create the parameter struct.
%
% Written by Julian Moosmann, 2017-06-02, last mod: 2017-06-03
%
% ADD_DATA_SET( restore_default_after_set, parameter_cell_array )

%% TODO: rename non-parameter variables (running index nn, default struct, etc)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    restore_default_after_set = 0;
end
if nargin < 2
    parameter_cell_array = evalin( 'caller', 'who');
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create struct from variables currently available in caller workspace
for nn = 1:numel( parameter_cell_array )
    name = parameter_cell_array{nn};
    
    if ~sum( strcmpi( name, {'PARAMETER_CELL', 'DEFAULT', 'ans', 'PARAMETER_STRUCT', 'raw', 'SUBSETS', 'DO_RECO', 'FIELDS_TO_PRINT'} ) )
        parameter_struct.(name) = evalin( 'caller', name );
    end
end

%% Initialize counter variable and define default
if ~isfield( parameter_struct, 'DATA_SET_NUM' )
    assignin( 'caller', 'DATA_SET_NUM', 0);
    parameter_struct.('DATA_SET_NUM') = 0;
    assignin( 'caller', 'DEFAULT', parameter_struct );   
elseif parameter_struct.DATA_SET_NUM == 0
    assignin( 'caller', 'DEFAULT', parameter_struct );   
end

%% Add parameter struct to struct array
assignin( 'caller', 'PARAMETER_STRUCT', parameter_struct);
evalin( 'caller', 'DATA_SET_NUM = DATA_SET_NUM + 1;' )
evalin( 'caller', 'PARAMETER_CELL{DATA_SET_NUM} = PARAMETER_STRUCT;' )
evalin( 'caller', 'PARAMETER_CELL{DATA_SET_NUM}.DATA_SET_NUM = DATA_SET_NUM;' )

%% Restore default parameters
if restore_default_after_set == 1
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
