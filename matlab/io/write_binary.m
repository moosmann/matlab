function write_binary(filename,im)
% write binary data
% filename: sting without ending

cl = class( im);
[s1,s2] = size(im);

fn = sprintf('%s_%ux%u.raw', filename, s1, s2);
fid = fopen(fn, 'w' );
fwrite(fid, im, cl );
fclose(fid );