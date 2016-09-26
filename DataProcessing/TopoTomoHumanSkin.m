% Folder: sl10p1hgx0p1vg_5minExT_SD0-1000_4steps
%% Read images
clear all
tic
if ~exist('rad','var')
    InputFolder = '/mnt/tomoraid-LSDF/tomo/TopoTomo_100927-SingleDistancePhaseRetrieval/radiography_human_skin/sl10p1hgx0p1vg_5minExT_SD0-1000_4steps/';
    filename = [InputFolder 'radio.tif'];
    rad = readmultitif(filename);
    [dimx dimy ~] = size(rad);
    roiX = ceil(dimx/4):floor(3*dimx/4);
    roiY = ceil(dimy/8):floor(7*dimy/8);
    rad  = rad(roiX,roiY,:);
    [dimx dimy NumIm] = size(rad);
    fprintf('Cropping data to %u x %u.\n',dimx,dimy)
    %fprintf('Read %u images of size %u x %u.\n',NumIm,dimx,dimy)
    % Filter hot pixels
    for ii = 1:NumIm
        rad(:,:,ii) = FilterHotPixel(rad(:,:,ii),0.05,0,[4 4]);
    end
end
%% Preprocessing parameters
NumDark = 5;
NumFlat = 3;
NumSam  = 1;
NumDist = (NumIm-NumDark)/(NumFlat+NumSam);
%% Dark field computation and correction
if ~exist('dark','var')
    dark = median(rad(:,:,1:NumDark),3);
    for ii = NumDark+(1:4*NumDist)
        rad(:,:,ii) = rad(:,:,ii)-dark;
    end
end
%% Flat field computation
if ~exist('flat','var')
    for ii = NumDist:-1:1
        flat(:,:,ii) = median(rad(:,:,NumDark+(ii-1)*(NumFlat+NumSam)+(1:NumFlat)),3);
    end
end
%% Flat field correction
if ~exist('int','var')
    for ii = NumDist:-1:1
        int(:,:,ii) = (rad(:,:,NumDark+(ii)*(NumFlat+NumSam)))./flat(:,:,ii);
    end
end
%% Phase retrieval
fprintf('Phase retrieval: ');tic
alphaCTF_alphaTIE = 2.5;
TIE_pCTF_CTF_PNLO = [1 0.1 0];
filterLowFreq = 60;
Padding_FactorAndValue = {2 'symmetric'};
for ii = NumDist:-1:1
    % intensity
    EnergyDistancePixelsize = [12 0.075+(ii-1)/3 1e-6];
    [out bf(:,:,ii)] = Reco(int(:,:,ii),alphaCTF_alphaTIE,EnergyDistancePixelsize,TIE_pCTF_CTF_PNLO,filterLowFreq,Padding_FactorAndValue);
    tie(:,:,ii)  = out.tieLO;
    pctf(:,:,ii) = out.ctfProjected;
    % flat field
    out = Reco(flat(:,:,ii),alphaCTF_alphaTIE,EnergyDistancePixelsize,TIE_pCTF_CTF_PNLO,filterLowFreq,Padding_FactorAndValue);
    flattie(:,:,ii)  = out.tieLO;
    flatpctf(:,:,ii) = out.ctfProjected; 
end
fprintf('Done in %gs.\n',toc)
%% Save images
OutputFolder = '/home/moosmann/TopoTomo_HumanSkin/edf/';
fprintf('Writing images to %s: ',OutputFolder);tic
filename = sprintf('%sdark.edf',OutputFolder);
edfwrite(filename,dark','float32')
bf = Binning(bf);
for ii = NumDist:-1:1
    filename = sprintf('%sint_z%04umm.edf',OutputFolder,round(1000*(0.075+(ii-1)/3)));
    edfwrite(filename,int(:,:,ii)','float32')
    filename = sprintf('%sphaseTIE_z%04umm.edf',OutputFolder,round(1000*(0.075+(ii-1)/3)));
    edfwrite(filename,tie(:,:,ii)','float32')
    filename = sprintf('%sphasePCTF_z%04umm.edf',OutputFolder,round(1000*(0.075+(ii-1)/3)));
    edfwrite(filename,pctf(:,:,ii)','float32')
    filename = sprintf('%sflat_z%04umm.edf',OutputFolder,round(1000*(0.075+(ii-1)/3)));
    edfwrite(filename,flat(:,:,ii)','float32')
    filename = sprintf('%sflatPhaseTIE_z%04umm.edf',OutputFolder,round(1000*(0.075+(ii-1)/3)));
    edfwrite(filename,flattie(:,:,ii)','float32')
    filename = sprintf('%sbinaryFilter_z%04umm.edf',OutputFolder,round(1000*(0.075+(ii-1)/3)));
    edfwrite(filename,bf(:,:,ii)','float32')
end
fprintf('Done in %gs\n',toc)