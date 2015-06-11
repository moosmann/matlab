% Cell of data sets
dataPath{1} = '/export/scratch1/moosmann/art/Batenburg/LowQualityDataSet/sino_HigherCounts';
dataPath{2} = '/export/scratch1/moosmann/art/Batenburg/LowQualityDataSet/sino_LowerCounts';
dataPath{3} = '/export/scratch1/moosmann/art/Batenburg/MediumQualityDataSet/sino';

% Loop over data sets
tic
fprintf(' Start Script\n')
for nn = 3%1:numel(dataPath)
    % Read sinogram filenames
    sinoName = FilenameCell(dataPath{nn});    
    % Loop over sinograms
    for mm = 1:numel(sinoName)
        % Read sino
        filename = sprintf('%s/%s',dataPath{nn},sinoName{mm});
        sino = imread(filename);
        % Loop over sino filter methods
        for ll = 4%1:3
            % Filter sino
            [sinof,filtStr] = FiltSino(sino,3,ll,MaskingTriangle(size(sino),210/148,[2 0]));
            % Write sino
            outputPath = sprintf('%s%s',dataPath{nn},filtStr);
            CheckAndMakePath(outputPath)
            filename = sprintf('%s/%s',outputPath,sinoName{mm});
            write32bitTIF(filename,sinof)
        end
    end
end
fprintf(' Finished in %g s = %.2g s\n',toc,toc/60)