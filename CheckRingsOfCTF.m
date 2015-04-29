function CheckRingsOfCTF(Energy,Distance,Pixelsize)

edp = [Energy Distance Pixelsize];
% Get dir struct and read image
imstruct = dir('*edf');
imname = imstruct(1).name;
int = pmedfread(imname)';
% Substrart mean: intensity contrast
g   = int - mean(int(:));
% Log of modulus of FT of intensity contrast
glaf  = laf(g);
% Projected CTF filter
f = SineFilter(size(int),edp,0.1);
% Plot.
itool(glaf)
%itool(f)
itool(glaf.*f)