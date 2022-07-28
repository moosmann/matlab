fn = '/beegfs/desy/user/moosmanj/drosophila/all_roi.txt';

fid = fopen(fn);
c = textscan(fid, '%s %s');


fclose(fid);

fn = '/beegfs/desy/user/moosmanj/drosophila/all_roi.txt';

fnt = '/beegfs/desy/user/moosmanj/drosophila/train_2.txt';
fnv = '/beegfs/desy/user/moosmanj/drosophila/val_2.txt';

fidt = fopen(fnt, 'w');
fidv = fopen(fnv, 'w');

for n = 1:size(c{1,1},1)
    
    if mod(n,5) == 0
        fprintf( fidv, '%s %s\n', c{1,1}{n}, c{1,2}{n} );
        fprintf( '\n %6u: train', n)
    else
        fprintf( fidt, '%s %s\n', c{1,1}{n}, c{1,2}{n} );
        fprintf( '\n %6u: valid', n)
    end
end

fclose(fidt);
fclose(fidv);