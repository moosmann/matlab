clear intd
%% Read data
InputFolder = '/mnt/tomoraid-LSDF/microscopy/NSRRC/20120422/output/sample13tomo1registered';
if ~exist('inty','var')
    inty = readstack(InputFolder,'Result*','tif',1,1,131,0.0);
end
if ~exist('dat','var')
    InputFolder = '/mnt/tomoraid-LSDF/microscopy/NSRRC/20120422/output/sample13tomo1';
    dat = readstack(InputFolder,'sample*','tiff',1,1,132,0.05);
    ref = dat(:,:,end);
    dat(:,:,end) = [];
    %ref = FilterHotPixel(double(imread([InputFolder 'sample 13-b2-90s-tomo1_refData.tiff'])),0.05);
    [dimx dimy NumIm] = size(dat);
end
% Read second flat.
%% ????????????????????????????????
filepath = '';
ref2 = double(imread(filepath));
refPos1 = 0;
%% ????????????????????????????????
refPos2 = 140;
%% Flat field correction
if ~exist('intd','var')
    for ii = NumIm:-1:1
        % Division
        %im = dat(:,:,ii)./ref;
        im = dat(:,:,ii)./squeeze(interp1(shiftdim(cat(3,refPos1,refPos2),1),shiftdim(cat(3,ref,ref2),2),ii,'linear','extrap'));
        %im = im - mean(im(:));
        intd(:,:,ii) = im;
        intdf(:,:,ii) = FilterDisk(im,3,1,0,0,2);
        % Substracton
%         im = dat(:,:,ii)-ref;
%         %im = im - mean(im(:));
%         ints(:,:,ii) = im;
%         intsf(:,:,ii) = FilterDisk(im,3,1,0,0,1);
    end
end
%% Image correlation
x = 150:350;
y = x;
for ii = (NumIm-1):-1:1
    out = ImageCorrelation(intdf(x,y,ii),intdf(x,y,ii+1),0,0);
    xshift(ii) = out.Xshift;
    yshift(ii) = out.Yshift;
end
figure('Name','Relative X shift VS projecton number'),plot(xshift)
figure('Name','Relative Y shift VS projecton number'),plot(yshift)
figure('Name','2D shift: relative X shift VS Relative Y shift'),plot(xshift,yshift,'x')

