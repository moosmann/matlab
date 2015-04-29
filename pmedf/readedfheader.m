% function hd=readedfheader(fid)

function hd=readedfheader(fid)

headerlength = 512;
closing = '}';
hd = fscanf(fid,'%c',headerlength);

while not(strcmp(hd(end-1),closing))
	hd = [hd fscanf(fid,'%c',headerlength)];
end
