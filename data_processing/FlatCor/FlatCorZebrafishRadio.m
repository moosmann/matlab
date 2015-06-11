ParentPath = '/mnt/tomoraid-LSDF/rci/MI1079-ESRF-ID19-July2011/inlineTomo/ZebraFish/ZebraFish_stage7d/ZebraFish_stage7d_12keV_radio_';
ParentPath = [ParentPath '/'];
DataFilePrefix = 'ZebraFish';
FlatStruct = dir([ParentPath 'flat*.edf']);
DataStruct = dir([ParentPath  DataFilePrefix '*.edf']);
HotPixThresh = 0.03;

for nn = numel(DataStruct):-1:1
    dat{nn} = FilterHotPixel(pmedfread([ParentPath DataStruct(nn).name]),HotPixThresh);
    flat{nn} = FilterHotPixel(pmedfread([ParentPath FlatStruct(nn+1).name]),HotPixThresh);
    
    int{nn} = dat{nn}./flat{nn};
end