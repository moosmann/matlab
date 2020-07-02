function assign_variable( name, value )
% Assign a default value of a variable or struct in the caller workspace if
% it does not exist.
%
% ARGUMENTS
% name : string, name of the variable or the struct to be checked
% value : any
%
% Written Julian Moosmann
%
% assign( name, value )

if ~isstruct( value )
    assignin( 'caller', 'assign_tmp', value )
    expr = sprintf( '%s = assign_tmp;', name  );
    evalin( 'caller', expr );
else
    fn = fieldnames( value );
    %for nn = 1:numel( fn )
        %fnn = fn{nn};
        fnn = fn{1};
        assignin( 'caller', 'assign_tmp', value.(fnn) )
        expr = sprintf( '%s.%s = assign_tmp;', name, fnn  );
        evalin( 'caller', expr );
        fprintf( '\nname: %s', name)
        fprintf( '\nname: %s', fnn)
    %end
end
