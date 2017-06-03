function ADD_DATA_SET( restore_default_after_set, cell_array_of_variable_names )
% Add (almost) all parameters, i.e. MATLAB variables given in the current
% workspace at the moment of calling this function to the cell struct
% 'par'.
%
% ARGUMENTS
% restore_default_after_set : bool. default: 0. restore parameter values
% back to default values after a data set is added to the loop.
% cell_array_of_variable_names : cell array of the names of the varibles to
% used to create the parameter struct.
%
% Written by Julian Moosmann, 2017-06-02, last mod: 2017-06-03
%
% ADD_DATA_SET( restore_default_after_set, cell_array_of_variable_names )

%% TODO: rename non-parameter variables (running index nn, default struct, etc)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    restore_default_after_set = 0;
end
if nargin < 2
    cell_array_of_variable_names = evalin( 'caller', 'who');
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create struct from variables currently available in caller workspace
for nn = 1:numel( cell_array_of_variable_names )
    name = cell_array_of_variable_names{nn};
    
    if ~sum( strcmpi( name, {'PARAMETER_STRUCT', 'default', 'ans', 'par', 'raw', 'nums', 'doreco', 'print_field'} ) )
        PARAMETER_STRUCT.(name) = evalin( 'caller', name );
    end
end

%% Initialize counter variable and define default
if ~isfield( PARAMETER_STRUCT, 'nn' )
    assignin( 'caller', 'nn', 0);
    PARAMETER_STRUCT.('nn') = 0;
    assignin( 'caller', 'default', PARAMETER_STRUCT );   
elseif PARAMETER_STRUCT.nn == 0
    assignin( 'caller', 'default', PARAMETER_STRUCT );   
end

%% Add parameter struct to struct array
assignin( 'caller', 'PARAMETER_STRUCT', PARAMETER_STRUCT);
evalin( 'caller', 'nn = nn + 1;' )
evalin( 'caller', 'par{nn} = PARAMETER_STRUCT;' )
evalin( 'caller', 'par{nn}.nn = nn;' )

%% Restore default parameters
if restore_default_after_set == 1
    % Get default struct from caller
    default = evalin( 'caller', 'default' );
    % Set parameter in caller back to defaults
    defn = fieldnames( default );
    for nn = 1:numel(defn)
        name = defn{nn};
        if ~strcmpi( name, 'nn' )
            assignin('caller', name, default.(name))
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
