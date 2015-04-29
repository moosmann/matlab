function out = testParsing

DataSetName = 'Xenos/stage28';
out = ParFunc('DataSet',DataSetName,'PixelRegion',{},'showFigure',0);
whos

function varargin = ParFunc(varargin)

%% PARSE INPUT ARGUMENT 'VARARGIN'
for nn = 1:2:numel(varargin)
    %eval([varargin{nn} ' = ' varargin{nn+1}]);
    %eval(varargin{nn}) = varargin{nn+1};
    assignin('caller',varargin{nn},varargin{nn+1});
end


