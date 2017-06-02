function parameter_struct = parse_workspace_variables_to_parameter_struct( cell_array_of_variable_names )
% Parse variables in caller workspace to a parameter struct. 
%
% ARGUMENTS
% cell_array_of_variable_names : cell array of the names of the varibles to
% used to create the parameter struct
%
% Written by Julian Moosmann, 2017-06-02
%
% parameter_struct = parse_workspace_variables_to_parameter_struct( cell_array_of_variable_names ) 

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    cell_array_of_variable_names = evalin( 'caller', 'who');
end
   
%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nn = 1:numel( cell_array_of_variable_names )
        name = cell_array_of_variable_names{nn};
        parameter_struct.(name) = evalin( 'caller', name );
    end
end