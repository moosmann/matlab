function var_value = assign_from_struct( parameter_struct, field_name, default_value, assign_struct )

%% Defaults
if nargin < 3
    default_value = [];
end
if nargin < 4
    assign_struct = 0;
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
    if assign_struct
        var_name = sprintf( '%s.%s', inputname(1), field_name );
        assignin( 'caller', 'assign_tmp', var_value )
        expr = sprintf( '%s = assign_tmp;', var_name  );
        evalin( 'caller', expr );
    else
        assignin( 'caller', field_name, var_value );
    end
end