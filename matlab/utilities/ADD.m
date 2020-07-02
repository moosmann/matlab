function ADD( how )
% Add parameters, i.e. MATLAB variables currently defined in the workspace
% at the caller function, to the global cell array of parameter structs
% 'PARAMETER_CELL'.
%
% ARGUMENTS
% how : string. default: 'data'. Defines how parameter sets should be added
%   'data' or '' : add current parameters set to the loop
%   'default' or 'd' : use current parameter set as default, does
%       not add a parameter set to the loop
%   'restore' or 'r' : add current parameter set to the loop and restore default
%       parameter afterwards. Thereby changes made after defining the
%       default parameter set are discarded.
%
% Written by Julian Moosmann, 2017-06-02, last mod: 2017-06-04
%
% ADD( how )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    how = 'data';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global DATA_SET_NUM
global DEFAULT
global PARAMETER_CELL

try
    % Retrieve names of variables currently defined in the caller workspace
    parameter_cell_array = evalin( 'caller', 'who');
    
    % Create struct from variables currently available in caller workspace
    for nn = 1:numel( parameter_cell_array )
        var_name = parameter_cell_array{nn};
        
        if ~sum( strcmpi( var_name, {'ans', 'raw', 'SUBSETS', 'RUN_RECO', 'PRINT_PARAMETERS'} ) )
            parameter_struct.(var_name) = evalin( 'caller', var_name );
        end
    end
    
    %% Initialize counter variable and default
    if ~isfield( parameter_struct, 'DATA_SET_NUM' )
        DATA_SET_NUM = 0;
        % workaround is needed to prevent counter increase when errors occur
        assignin( 'caller', 'DATA_SET_NUM', DATA_SET_NUM )
        DEFAULT = parameter_struct;
        PARAMETER_CELL = [];
    end
    
    switch lower( how )
        
        case {'', 'data'}
            %% Add parameter set to cell array of parameter structs
            
            DATA_SET_NUM = DATA_SET_NUM + 1;
            parameter_struct.DATA_SET_NUM = DATA_SET_NUM;
            PARAMETER_CELL{DATA_SET_NUM} = parameter_struct;
            
        case {'d', 'default'}
            %% Set parameter set as default
            DEFAULT = parameter_struct;
            
        case {'r', 'restore'}
            %% Add parameter set AND Restore default parameters
            
            DATA_SET_NUM = DATA_SET_NUM + 1;
            parameter_struct.DATA_SET_NUM = DATA_SET_NUM;
            PARAMETER_CELL{DATA_SET_NUM} = parameter_struct;
            
            % Set parameter in caller back to defaults
            fn = fieldnames( DEFAULT );
            for nn = 1:numel( fn )
                var_name = fn{nn};
                assignin('caller', var_name, DEFAULT.(var_name))
            end
    end
catch exception
    clear global
    rethrow( exception )
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
