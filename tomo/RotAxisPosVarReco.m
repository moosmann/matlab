function [vol RotPos] = RotAxisPosVarReco(sino,RotAxisPosGuess,PosHalfRange_PosInc,HorPixHalfRange,scanRangeInDeg)

if nargin < 2
    RotAxisPosGuess = round(size(sino,2)/2);
end
if nargin < 3
    PosHalfRange_PosInc = [100 10];
end
if nargin < 4
    HorPixHalfRange = 0;
end
if nargin < 5
    scanRangeInDeg = 180;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotation axis position scan range
PosHalfRange = PosHalfRange_PosInc(1);
if numel(PosHalfRange_PosInc) == 2
    PosInc = PosHalfRange_PosInc(2);
else
    PosInc = 1;
end
RotPos = RotAxisPosGuess + (-PosHalfRange:PosInc:PosHalfRange);
% horizontal detector range
if HorPixHalfRange == 0
    HorPixHalfRange = floor( size(sino,2)/2 -PosHalfRange_PosInc(1) -1);
end    
HorPixRange = -HorPixHalfRange:HorPixHalfRange;
% angular range
NumProjTomo = size(sino,1);
angles = scanRangeInDeg/180*pi/(NumProjTomo)*((1:NumProjTomo)-1);

%% loop
for nn = numel(RotPos):-1:1
    fprintf('%6u',nn) 
end
fprintf('\n')
for nn = numel(RotPos):-1:1
    fprintf('%6u',RotPos(nn))
    vol(:,:,nn) = astra_make_reco(SubtractMean(sino(:,RotPos(nn) + HorPixRange)),angles,'FBP_CUDA',1);
end
fprintf('\n')
x0 = ceil(size(vol,1)/sqrt(2));
x = floor((size(vol,1)-x0)/2)+(1:x0);
nimplay(vol(x,x,:))
