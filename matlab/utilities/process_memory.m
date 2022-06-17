function m = process_memory( process_id )
% Returns memory consumption in kilobyte in for process ID.
    s = sprintf( 'pmap %u | tail -n 1', process_id);
    
    [~,t] = unix( s );
    c = textscan( t, '%s %uK');
    m = c{2};
    