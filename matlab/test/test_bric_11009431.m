f = dir('/asap3/petra3/gpfs/p07/2020/data/11009431/raw/bric*');

for nn = 1:numel(f)
    
    ff = [f(nn).folder filesep f(nn).name filesep 'tiff*/*tif'];
    d = dir( ff);
    a = [d(:).bytes];
    num_kaput = 0;
    if ~isempty( a)
        m = a == a(1);
        num_kaput = sum(m==0);        
        c = {d(~m).name};
        if num_kaput
%             keyboard
        end
    end    
    fprintf( '\n%2u: %s %u %u', nn, f(nn).name, numel(d), num_kaput )
    
%     for mm = 1:numel(c)
%         fprintf( '\n  %s', c{mm})
%     end
    ind = 1:numel(d);
    fprintf( '\n  ')
    indk = ind(~m);
    fprintf(' %u', indk)
    fprintf( '\n')
    fprintf(' %u', indk(2:end) -indk(1:end-1) )
    fprintf( '\n')
end