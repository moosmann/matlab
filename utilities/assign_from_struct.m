function var_value = assign_from_struct( parameter_struct, field_name, default_value )

%% Defaults
if nargin < 3
    default_value = [];
end

%% Main
if isfield( parameter_struct, field_name )    
    var_value = parameter_struct.( field_name );
    if isempty( var_value )
        var_value = default_value;
    end
else
    var_value = default_value;
end

if nargout == 0
    assignin( 'caller', field_name, var_value );
end