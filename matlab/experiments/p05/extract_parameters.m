fid =fopen('p05_reco.m','r');

while 1
    [c, pos] = textscan( fid, '%s', 1, 'Delimiter', {'\n', '\r'} );
    s = c{1}{1};
    if strcmpi( s, '%% PARAMETERS / SETTINGS %%' )
        tag = 1;
    end
    
    if strcmpi( s, '%% END OF PARAMETERS / SETTINGS %%' )
        c_break = c;
        break
    end
    
    if tag
        
    end
end

c{1}{1}

fclose(fid);


