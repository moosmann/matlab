% Experiment at TopoTomo@ANKA, 2012-08-27, Radiography of Cells

InputPath = '/mnt/tomoraid-LSDF/UserData/RadiographyOfCells_20120827/';
%InputPath = '/mnt/LSDF/UserData/RadiographyOfCells_20120827';
HotPixThreshold = 0.02;
printInfo = true;

%% Read images
% dark0 seems to be less noise than dark1. Mean value and std are almost
% identical. Min/Max variation differ by 30 %
if ~exist('dark0','var')
    % First day, evening/night
    dark0 = FilterHotPixel(double(imread([InputPath 'sample_dark.tif'])),HotPixThreshold,printInfo);
    flat0 = FilterHotPixel(double(imread([InputPath 'sample_flat.tif'])),HotPixThreshold,printInfo)-dark0;
    radio0 = FilterHotPixel(double(imread([InputPath 'sample_radio.tif'])),HotPixThreshold,printInfo)-dark0;
    % Second day, morning
    dark1 = FilterHotPixel(double(imread([InputPath 'L8_dark_2minutes.tif'])),HotPixThreshold,printInfo);
    flat1 = FilterHotPixel(double(imread([InputPath 'L8_flat_2minutes.tif'])),HotPixThreshold,printInfo)-dark1;
    flat2 = FilterHotPixel(double(imread([InputPath 'L8_flat_2minutes_2.tif'])),HotPixThreshold,printInfo)-dark1;
    radio1 = FilterHotPixel(double(imread([InputPath 'L8_radio_2minutes.tif'])),HotPixThreshold,printInfo)-dark1;
    radio2 = FilterHotPixel(double(imread([InputPath 'L8_radio_2minutes_2.tif'])),HotPixThreshold,printInfo)-dark1;
    
    % test = FilterHotPixel(double(imread('test.tif')),HotPixThreshold,printInfo);
    % test2 = FilterHotPixel(double(imread('test2.tif')),HotPixThreshold,printInfo);
    % test3 = FilterHotPixel(double(imread('test3.tif')),HotPixThreshold,printInfo);
    
    %% Standard flat field correction
    int0 = radio0./flat0;
    int1 = radio1./flat1;
    int2 = radio2./flat2;
    
    %% Print domains
    domain(dark0),domain(flat0),domain(radio0)
    domain(dark1),domain(flat1),domain(radio1)
    domain(flat2),domain(radio2)
end


%% Crop ROI for image correlation
% Image 0
rpos0 = [1760         975         160         134];
fpos0 = [2158         842         160         134];
xf0 = fpos0(2) + (1:fpos0(4));
yf0 = fpos0(1) + (1:fpos0(3));
xr0 = rpos0(2) + (1:rpos0(4));
yr0 = rpos0(1) + (1:rpos0(3));
out0 = ImageCorrelation(radio0(xr0,yr0),flat0(xf0,yf0),1);
d0 = abs(rpos0-fpos0);
flat0c  = flat0(1:end-d0(2),d0(1)+1:end);
radio0c = radio0(d0(2)+1:end,1:end-d0(1));
% Image 1
rpos1 = [1187         218         307         191];
%fpos1 = rpos1 + [d0(1) -d0(2) 0 0];
fpos1 = [1595          85         307         191];
xf1 = fpos1(2) + (1:fpos1(4));
yf1 = fpos1(1) + (1:fpos1(3));
xr1 = rpos1(2) + (1:rpos1(4));
yr1 = rpos1(1) + (1:rpos1(3));
out1 = ImageCorrelation(radio1(xr1,yr1),flat1(xf1,yf1),1);
d1 = abs(rpos1-fpos1);
flat1c  = flat1(1:end-d1(2),d1(1)+1:end);
radio1c = radio1(d1(2)+1:end,1:end-d1(1));
% Image 2
rpos2 = [1209         182         328         217];
fpos2 = [1196         184         328         217];
xf2 = fpos2(2) + (1:fpos2(4));
yf2 = fpos2(1) + (1:fpos2(3));
xr2 = rpos2(2) + (1:rpos2(4));
yr2 = rpos2(1) + (1:rpos2(3));
out2 = ImageCorrelation(radio2(xr2,yr2),flat2(xf2,yf2),1);
d2 = abs(rpos2-fpos2);
flat2c  = flat2(1:end-d2(2),d2(1)+1:end);
radio2c = radio2(d2(2)+1:end,1:end-d2(1));

%% Shifted flat field correction
% Image 0
ipos0 = [2009         922         859         740];
x = ipos0(2) + (1:ipos0(4));
y = ipos0(1) + (1:ipos0(3));
im0 = radio0c./flat0c;
im0c = im0(x,y);
% Image 1
ipos0 = [833        1182        1344         879];
x = ipos0(2) + (1:ipos0(4));
y = ipos0(1) + (1:ipos0(3));
im1 = radio1c./flat1c;
im1c = im1(x,y);

%out = Reco(im0c,2.5,[15 0.187 0.36e-6],[1 0.1]);
%% Phase retrieval
if true
    nf = 20;
    for nn = nf:-1:1
        out = RecoGPU(im1c,0.1+2*(nn-1)/(nf-1),[15 0.187 0.36e-6],0,0,0,0.1);
        %domain(ctf(:,:,nn))
        ctf1(:,:,nn) = RemoveLowFreq(out.ctfProjected);
    end
    for nn = nf:-1:1
        out = RecoGPU(im0c,0.1+2*(nn-1)/(nf-1),[15 0.187 0.36e-6],0,0,0,0.1);
        %domain(ctf(:,:,nn))
        ctf0(:,:,nn) = RemoveLowFreq(out.ctfProjected);
    end
    % out1 = ImageCorrelation(radio1(x,y),flat1(x,y),1,1);
    % out2 = ImageCorrelation(radio2(x,y),flat2(x,y),1,1);
end