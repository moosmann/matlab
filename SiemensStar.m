function star = SiemensStar(resolution,spokes)
% Siemens star like pattern. Input paramters: resolution, (number of)
% spokes.
    
    if nargin<1,resolution=[512,512];end;
    if nargin<2,spokes=16;end;
    
if size(resolution,2)==1,
    dimx = resolution;
    dimy = resolution;
elseif size(resolution,2)==2,
    dimx = resolution(1);
    dimy = resolution(2);
else
    fprintf(1,'Error resolution is not a 1x2 matrix\n');
end;

star   = ones(dimx+1,dimy+1);
disc   = zeros(dimx+1,dimy+1);
center = disc;
%[x,y]  = meshgrid(-1/2:1/dimx:1/2-1/dimx,-1/2:1/dimy:1/2-1/dimy);
[x,y]  = meshgrid(-1/2:1/dimx:1/2,-1/2:1/dimy:1/2);
theta  = mod(pi/2+atan(x./y),4*pi/spokes);
disc((sqrt(x.^2+y.^2)>0.4))=1;
%center((sqrt(x.^2+y.^2)<spokes/4/dimx))=1; 
center((sqrt(x.^2+y.^2)<spokes/4/512))=1; 

if mod(spokes,4)>0,
   theta = cat(1,theta(1:dimx/2,:),theta(dimx/2+1:end,end:-1:1));
end;

star((disc|((spokes*theta>=3*pi))|(spokes*theta<=1*pi))) = 0;
star(center>0) = 1;
star = star(1:dimx,1:dimy);
