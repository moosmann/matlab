function data = renormalize(name,nop);
    
dimx = 1024;
dimy = 1024;


data = zeros(dimy,dimx,nop);

first = 1000;last =first+ nop;
for ii=first:last
    
%get filename string
namestr = sprintf('%s%4u.edf',name,ii);

%read data and header
data(:,:,ii) = pmedfread(namestr);

end;

