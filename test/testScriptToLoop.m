scan_path = 'nix';
para2 = 0;

if exist( 'external_parameter' ,'var')
    fprintf('\nexternal paramter exists\n')
    
    fields = fieldnames( external_parameter );
    
    for nn = 1:numel(fields)
        var_name = fields{nn};
        var_val = getfield(external_parameter, var_name );
        assignin('base', var_name, var_val )
    end
end



fprintf('\nscan_path : %s', scan_path )
fprintf('\npara2 : %g\n\n', para2 )