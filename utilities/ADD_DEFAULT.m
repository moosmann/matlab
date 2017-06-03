function ADD_DEFAULT( parameter_cell_array )
% Add all parameters, i.e. MATLAB variables given in the current
% workspace at the moment of calling this function to the cell struct
% 'DEFAULT'.
%
% ARGUMENTS
% parameter_cell_array : cell array of the names of the varibles to
% used to create the parameter struct.
%
% Written by Julian Moosmann, 2017-06-02, last mod: 2017-06-03
%
% ADD_DEFAULT( parameter_cell_array )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    parameter_cell_array = evalin( 'caller', 'who');
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create struct from variables currently available in caller workspace
for nn = 1:numel( parameter_cell_array )
    name = parameter_cell_array{nn};
    
    if ~sum( strcmpi( name, {'SUBSETS', 'RUN_RECO', 'PRINT_PARAMETERS'} ) )
        parameter_struct.(name) = evalin( 'caller', name );
    end
end

%% Set default
assignin( 'caller', 'DATA_SET_NUM', 0);
parameter_struct.('DATA_SET_NUM') = 0;
assignin( 'caller', 'DEFAULT', parameter_struct );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
