function CheckRotAxis(im0,im180,xWidth,yRange)

if nargin  < 3
    xWidth = 25;
end
if nargin < 4
    yRange = 400:1200;
end
% % get first and last projecton from stack
% im  = squeeze(stack(:,:,1));
% im2 = squeeze(stack(:,:,end));

% Loop along vertical direciton
nmax = floor(size(im0,1)/xWidth)-1;
rr(nmax+1) = 0;
fprintf('\nImages size: %u x %u:\n',size(im0))
fprintf('Rotation axis position along horizontal sections of width %u:\n',xWidth)
for nn = 0:nmax
    x = nn*xWidth +(1:xWidth);
    out = ImageCorrelation(im0(x,yRange),fliplr(im180(x,yRange)),0,0);
    RotAxisPos = out.VerticalRotationAxisPosition;
    rr(nn+1) = RotAxisPos;
    fprintf('%4u: %6.1f  ',x(1),RotAxisPos);
    if ~mod((nn+1),10)
        fprintf('\n');
    end
end
fprintf('\n')
figure('Name',sprintf('RotAxisPos VS z. median(RotAxisPos)=%.1f, xWidth=%u, yRange(1)=%u:%u',median(rr),xWidth,yRange(1),yRange(end)))
plot(rr)