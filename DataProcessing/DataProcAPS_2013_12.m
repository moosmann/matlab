% Script to process data acquired at 2-BM-B at APS in November 2013,
% intended to test Neural-Network reconstructions. The training data sets
% are however shifted and rotated w.r.t. each other.
%
% First images of first scan are dark.
%
% Written by Julian Moosmann, last version 2014-01-30

% sample    proj_2
% rotaxis   1016

clear all
%% Parameter
printInfo = 0;
offset = 1;
readPhase(1) = 1;
readInt(1) = 1;
normInt(1) = 1;
filtSino(1) = 1;
%RotAxisPos = [1016 1016 1030 1016]; manually
RotAxisPos = [1013 1013 1026.5 1028];
% folder
ParentPath = '/mnt/LSDF/tomo/APS/2013-12_Xenopus';
DataSets = FilenameCell(sprintf('%s/data/',ParentPath));
if normInt
    phaseStr = 'phaseFluc';
else
    phaseStr = 'phase';
end
if filtSino
    appStr = '_filtSino';
else
    appStr = '';
end

%% Loop over data
for nn = 1:numel(DataSets)
    SubSets = FilenameCell(sprintf('%s/data/%s/proj*',ParentPath,DataSets{nn}));
    for ll = 1:numel(SubSets)
        % input folder
        InputPath = sprintf('%s/data/%s/%s',ParentPath,DataSets{nn},SubSets{ll});
        fprintf('\nProcessing data set: %s',InputPath)
        % create output folder
        OutputPathInt = sprintf('%s/int/%s/%s',ParentPath,DataSets{nn},SubSets{ll});
        CheckAndMakePath(OutputPathInt)
        OutputPathPhase = sprintf('%s/%s/%s/%s',ParentPath,phaseStr,DataSets{nn},SubSets{ll});
        CheckAndMakePath(OutputPathPhase)
        OutputPathSinoPhase = sprintf('%s/sino%s/%s/%s/%s',ParentPath,appStr,phaseStr,DataSets{nn},SubSets{ll});
        CheckAndMakePath(OutputPathSinoPhase)
        OutputPathVolPhase = sprintf('%s/vol%s/%s/%s/%s',ParentPath,appStr,phaseStr,DataSets{nn},SubSets{ll});
        CheckAndMakePath(OutputPathVolPhase)
        % filenames
        darkNames = FilenameCell(sprintf('%s/*Dark*',InputPath));
        flatNames = FilenameCell(sprintf('%s/*Flat*',InputPath));
        projNames = FilenameCell(sprintf('%s/%s_*',InputPath,SubSets{ll}));
        numDark = numel(darkNames);
        numFlat = numel(flatNames);
        
        numProj = numel(projNames);
        % get images size
        filename = sprintf('%s/%s',InputPath,projNames{2});
        iminfo = imfinfo(filename,'tif');
        dim1 = iminfo.Height;
        dim2 = iminfo.Width;
        % create phase retrieval filter
        phaseFilter = PhaseFilter('tie',[dim1 dim2],[22.5 0.670 1.1e-6],2.5);
        %% Processing
        if ~readInt
            % darks
            fprintf('\n Processing darks')
            stack = zeros(dim1,dim2,numDark-offset);
            parfor kk = 1:numel(darkNames)-offset
                filename = sprintf('%s/%s',InputPath,darkNames{kk+offset});
                stack(:,:,kk) = FilterPixel(single(imread(filename)),[0.05 0.01],printInfo);
            end
            darkMean = median(stack,3);
            % flats
            fprintf('\n Processing flats')
            stack = zeros(dim1,dim2,numFlat-offset);
            parfor kk = 1:numel(flatNames)-offset
                filename = sprintf('%s/%s',InputPath,flatNames{kk+offset});
                stack(:,:,kk) = FilterPixel(single(imread(filename)),[0.03 0.005],printInfo);
            end
            flatMean = median(stack,3);
            % some pixels in upper right corner are still zeros
            flatMean(flatMean == 0) = mean(flatMean(:));
            % mask microscope aperture
            %mask = MaskingAperture(darkMean,flatMean,10);
            %mask = floor(mask);
            mask = MaskingDisc(ones(size(darkMean)),0.99,0);
            darkMean = darkMean.*mask;
            flatMean = flatMean - darkMean;
            % invert flat
            %flatMean = FilterPixel(1./flatMean,[0.01,0.005],1);
            flatMean = 1./flatMean;
        end
        % projections
        fprintf('\n Processing projections\n')
        stack = zeros(dim1,dim2,numProj);
        if readPhase
            for kk = 1:numel(projNames)
                PrintNum(kk);
                % read phase map
                filename = sprintf('%s/phase_%04u.tif',OutputPathPhase,kk);
                im = imread(filename);
                % append image to stack
                stack(:,:,kk) = im;
            end
        else
            for kk = 1:numel(projNames)
                PrintNum(kk);
                if readInt
                    % read int from disk
                    filename = sprintf('%s/int_%04u.tif',OutputPathInt,kk);
                    im = imread(filename);
                else
                    % read proj from disc and filter pixel
                    filename = sprintf('%s/%s',InputPath,projNames{kk});
                    im = FilterPixel(single(imread(filename)),[0.03 0.005],printInfo);
                    %% flat/dark field correction
                    im = (im - darkMean).*flatMean;
                    % aperture masking
                    im = im.*mask;
                    im = im + sum(im(:))/sum(mask(:))*(1-mask);
                    % write intensities
                    filename = sprintf('%s/int_%04u',OutputPathInt,kk);
                    WriteImage(filename,im,'tif')
                    
                end
                %% phase retrieval
                %normalize projection
                if normInt
                    im = im/mean(im(:));
                end
                im = real( ifft2( phaseFilter.* fft2( im - 1 ) ) );
                % write phase maps
                filename = sprintf('%s/phase_%04u',OutputPathPhase,kk);
                WriteImage(filename,im,'tif')
                % append image to stack
                stack(:,:,kk) = im;
            end
        end
        %% tomographic reconstruction
        fprintf('\n Tomographic reconstruction\n')
        angles = 1*pi/numProj * ((1:numProj)-1);
        for kk = dim1:-1:1
            PrintNum(kk);
            % Filter sinogramm
            im = squeeze(stack(kk,:,:));
            if filtSino
                % last arguments: 1 means vertical, 2 means horizontal
                im = FilterSino(im,2,1);
            end
            % Create and write sino
            sino = SubtractMean(RotAxisSymmetricCropping(im',RotAxisPos(nn+(ll-1)) ));
            filename = sprintf('%s/sino_%04u',OutputPathSinoPhase,kk);
            WriteImage(filename,sino,'tif');
            % Reconstruct tomographic slice
            vol = astra_make_reco(sino,angles,'FBP_CUDA',1);
            % Write tomo slice
            filename = sprintf('%s/slice_%04u',OutputPathVolPhase,kk);
            WriteImage(filename,vol,'tif');
        end
        clear stack;
    end
end
fprintf('\nFinished \n\n')