function Embryo24SlitsClosedDistScan(PixelRegion)
   
% Parameters.
NumDarks = 5;
NumFlats = 3;
NumData = 1;

% PROGRAMME.    
% Path to data.
folder = '/mnt/tomoraid3/tomo/TopoTomo_SinglePhaseRetrieval_Jan2011/';
%scan = 'radiography/glue_wax_distances/radio_0p09vx0h_sl1.tif';
% scan = 'radiography_embryo/stage24_slits1_0p1vx0p1h_ET10min_SD0-1000/radio.tif';
scan = 'radiography_embryo/stage24_slits1_2vx2h_ET1p5min_SD0-1000/radio.tif';
% Read multi tif file.
stack = readmultitif([folder scan])/4.0;
if nargin>0
    PixelRegion = [351 1872; 101 3908];
    stack = stack(PixelRegion(1,1):PixelRegion(1,2),PixelRegion(2,1):PixelRegion(2,2),:);
    fprintf('Size of cropped stack: %u %u %u\n',size(stack))
end
%assignin('base',['stack'],stack);
% Assign dark, flat and data images.
[dimx dimy NumTotal] = size(stack);
NumPos  = (NumTotal-NumDarks)/(1+NumFlats);
darks   = stack(:,:,1:NumDarks);
assignin('base','darks',darks);
%dark   = filter_im(median(stack(:,:,1:NumDarks),3));
dark   = median(stack(:,:,1:NumDarks),3);

flat = zeros(dimx,dimy,NumPos);
for n=1:NumPos,
    flat(:,:,n) = median(stack(:,:,NumDarks+(NumFlats+NumData)*(n-1)+(1:NumFlats)),3);
end;

% Make variables available in user's workspace.
assignin('base','flat',flat);
assignin('base','flats1',stack(:,:,NumDarks+(NumFlats+NumData)*(1-1)+(1:NumFlats)));
assignin('base','flats2',stack(:,:,NumDarks+(NumFlats+NumData)*(2-1)+(1:NumFlats)));
assignin('base',['flats' num2str(n)],stack(:,:,NumDarks+(NumFlats+NumData)*(n-1)+(1:NumFlats)));
data   = stack(:,:,NumDarks+1+NumFlats:1+NumFlats:end);
assignin('base','data',data);
int    = (data - repmat(dark,[1 1 NumPos]))./(flat - repmat(dark,[1 1 NumPos]));
assignin('base','int',int);
int1    = (data)./(flat);
assignin('base','int1',int1);
assignin('base','dark',dark);








