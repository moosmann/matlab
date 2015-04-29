function ScriptExpTopotomo()
   
if nargin<1
    PixelRegion = 0;
end;   
 
% Parameters.
NumDarks = 5;
NumFlats = 3;
PixelRegion = {[101 2172] [101 3908]};

% PROGRAMME.    
% Path to data.
folder = '/mnt/tomoraid3/tomo/TopoTomo_SinglePhaseRetrieval_Jan2011/';
%scan   = 'radiography/glue_wax_distances/radio_0p09vx0h_sl1.tif';
scan   = 'radiography_embryo/stage24_slits1_0p1vx0p1h_ET10min_SD0-1000/radio.tif';
% Read multi tif file.
if size(PixelRegion,2)==1
    stack = readmultitif([folder scan])/4.0;
else
    stack = readmultitif([folder scan],PixelRegion)/4.0;
end;
assignin('base',['stack'],stack);
%stack  = readmultitif([folder scan],PixelRegion)/4;
% Assign dark, flat and data images.
[~ , ~, N] = size(stack);
N      = (N-NumDarks)/2;
dark   = stack(:,:,1:NumDarks);
assignin('base',['darks'],dark);
dark   = sum(dark,3)/NumDarks;
%flat   = stack(:,:,NumDarks+1:2:end-1)-repmat(dark,[1 1 N]);
%data   = stack(:,:,NumDarks+2:2:end)-repmat(dark,[1 1 N]);
flat   = stack(:,:,NumDarks+1:2:end-1);
data   = stack(:,:,NumDarks+2:2:end);
int    = (data)./(flat);

% Make variables available in user's workspace.
assignin('base',['dark'],dark);
assignin('base',['flat'],flat);
assignin('base',['data'],data);
assignin('base',['int'],int);
