function Embryo24SlitsClosedDistScanCell(PixelRegion)
   
if nargin<1
    PixelRegion = {[351 1872] [101 3908]};
    %PixelRegion = {[601 1672] [1001 2908]};
end;   
 
% Parameters.
NumDarks = 5;
NumFlats = 3;

% PROGRAMME.    
% Path to data.
folder = '/mnt/tomoraid3/tomo/TopoTomo_SinglePhaseRetrieval_Jan2011/';
%scan   = 'radiography/glue_wax_distances/radio_0p09vx0h_sl1.tif';
scan   = 'radiography_embryo/stage24_slits1_0p1vx0p1h_ET10min_SD0-1000/radio.tif';
% Read multi tif file.
if size(PixelRegion,2)==1
    stack = readMultiTif([folder scan],1/4);
else
    stack = readMultiTif([folder scan],PixelRegion,1/4.0);
end;
%assignin('base',['stack'],stack);
% Assign dark, flat and data images.
N = numel(stack);
N      = (N-NumDarks)/(1+NumFlats);
dark   = median(cat(3,stack{1:NumDarks}),3);
flat = cell(1,N);
for n=1:N,
    flat{n} = median(cat(3,stack{NumDarks+4*(n-1)+(1:NumFlats)}),3);
end;

% Make variables available in user's workspace.
assignin('base','flat',flat);
data   = stack(:,:,NumDarks+1+NumFlats:1+NumFlats:end);
assignin('base','data',data);
int    = (data - repmat(dark,[1 1 N]))./(flat);
assignin('base','int',int);
int1    = (data)./(flat);
assignin('base','int1',int1);
assignin('base','dark',dark);








