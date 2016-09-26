

ParentPath  = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/data/';
FolderStruct = dir([ParentPath '*mm*']);
fprintf('\nDATA SETS:\n\n')
for nn = 1:numel(FolderStruct)
    fprintf('%2u: %s\n',nn,FolderStruct(nn).name)
end
DataSetNum = input('\nType in number of the data set to reconstruct: ');
fprintf('\nDATA SET WHICH WILL BE RECONSTRUCTED:\n\n%2u: %s\n',DataSetNum,FolderStruct(DataSetNum).name)