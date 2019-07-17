function assign_default( name, default_value, assign_if_empty )
% Assign a default value of a variable or struct in the caller workspace if
% it does not exist.
%
% ARGUMENTS
% name : string, name of the variable or the struct to be checked
% default_value : any
% assign_if_empty : bool
%
% Written Julian Moosmann
%
% assign_default( name, default_value, assign_if_empty ) 

%% Default arguments
if nargin < 2
    default_value = [];
end
if nargin < 3
    assign_if_empty = 0;
end

%% Variable or struct with fields
r = regexp( name, '\.' );
if isempty( r )
    expr = sprintf( 'exist( ''%s'', ''var'' )', name );
    var_exists = evalin( 'caller', expr );
    name_exists = var_exists;
else
    %% Go down struct
    var_name = name(1:r(1)-1);
    expr = sprintf( 'exist( ''%s'', ''var'' )', var_name );
    var_exists = evalin( 'caller', expr );
    
    if ~var_exists
        name_exists = 0;
    else
        num_char = numel( name );
        num_subs = numel( r );
        % Check if field of struct exist
        if var_exists
            for nn = 1:num_subs
                sub_struct_name = name(1:r(nn)-1);
                if nn == num_subs
                    field_name = name( r(nn) + 1:num_char );
                else
                    field_name = name( r(nn) + 1 : r(nn+1) -1 );
                end
                expr = sprintf( 'isfield( %s, ''%s'')', sub_struct_name, field_name );
                name_exists = evalin( 'caller', expr );
                if ~name_exists
                    break
                end
            end
        end
    end
end

%% If name exists, check if empty
name_is_empty = 0;
if name_exists && assign_if_empty
    expr = sprintf( 'isempty( ''%s'')', name );
    name_is_empty = ~evalin( 'caller', expr );
end

%% Assignment via temp variable
if ~name_exists || (name_exists && assign_if_empty && name_is_empty)
    assignin( 'caller', 'assign_tmp', default_value )
    expr = sprintf( '%s = assign_tmp;', name  );
    evalin( 'caller', expr );
end
