function [fbp, art] = testAstraSaveImages(fbp, artStack, error, ParentPath,recType,showFigs)
% Function for 'testAstra' to create and save figures and images of
% reconstruction and their differences.
%
% Written by Julian Moosmann, last version 2013-10-28


if nargin < 5
    showFigs = 1;
end

%% Defaults and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MakePath(ParentPath);
% Image output format
[NumPix, ~, NumIter] = size(artStack); 
% roi range
dx = ceil(NumPix*(2-sqrt(2))/4*1.1);
x = dx:NumPix-dx;

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimum of error to FBP full
%[~,fbpErrMinPos] = min(error.full.fbp.full);
[~,fbpErrMinPos] = min(error.fbp.full);
if fbpErrMinPos < NumIter
    fprintf(' Position of minimum of error ART-FBP_full: iteration %u of %u\n',fbpErrMinPos,NumIter)
    art.opt = artStack(:,:,fbpErrMinPos);
elseif fbpErrMinPos == NumIter
    fprintf(' Position of minimum of error ART-FBP_full: last iteration %u\n',NumIter)
else
    fprintf(' ERROR: Positon of minimum %u is larger that number of iterations %u.\n',fbpErrMinPos,NumIter)
end
if NumIter ~= numel(error.res)
    fprintf(' \nn!! Size of variables (number of iterations) does not match!!\n\n')
end
art.last = artStack(:,:,NumIter);
fprintf(' Output path: %s\n',ParentPath)
%% Errors
switch showFigs
    case 1
        VisToggle = 'on';
    case 0
        VisToggle = 'off';
end

filename = sprintf('error__%s_resdiual',recType);
figure('Visible',VisToggle,'Name',sprintf('Error: %s',filename)); 
plot(error.res,'.')
saveas(gcf,sprintf('%s%s.eps',ParentPath,filename),'epsc2')

filename = sprintf('error__%s-fbp_full',recType);
figure('Visible',VisToggle,'Name',sprintf('Error: %s',filename)); 
plot(error.fbp.full,'.')
saveas(gcf,sprintf('%s%s.eps',ParentPath,filename),'epsc2')

filename = sprintf('error__%s-fbp_red',recType);
figure('Visible',VisToggle,'Name',sprintf('Error: %s',filename)); 
plot(error.fbp.red,'.')
saveas(gcf,sprintf('%s%s.eps',ParentPath,filename),'epsc2')


%% Full images
prefix = sprintf('%sreco_full__',ParentPath);
postfix = '';
imFormat = 'tif';
% Recos
WriteImage( sprintf('%sfbp_full%s',prefix,postfix), fbp.full, imFormat)
WriteImage( sprintf('%sfbp_red%s',prefix,postfix), fbp.red, imFormat)
WriteImage( sprintf('%sfbp_blur%s',prefix,postfix), fbp.blur, imFormat)
WriteImage( sprintf('%s%s_last__iter%04u%s',prefix,recType,NumIter,postfix), art.last, imFormat)
% Differences
WriteImage(sprintf('%sdiff__fbp_red-fbp_full%s',prefix,postfix), fbp.red-fbp.full, imFormat)
WriteImage( sprintf('%sdiff__fbp_blur-fbp_full%s',prefix,postfix), fbp.blur-fbp.full, imFormat)
WriteImage( sprintf('%sdiff__%s_last-fbp_full__iter%04u%s',prefix,recType,NumIter,postfix), art.last-fbp.full, imFormat)
WriteImage( sprintf('%sdiff__%s_last-fbp_blur__iter%04u%s',prefix,recType,NumIter,postfix), art.last-fbp.blur, imFormat)
if fbpErrMinPos < NumIter
    WriteImage( sprintf('%s%s_opt__iter%04u%s',prefix,recType,fbpErrMinPos,postfix), art.opt, imFormat)
    WriteImage( sprintf('%sdiff__%s_opt-fbp_full__iter%04u%s',prefix,recType,fbpErrMinPos,postfix), art.opt-fbp.full, imFormat)
    WriteImage( sprintf('%sdiff__%s_opt-fbp_blur__iter%04u%s',prefix,recType,fbpErrMinPos,postfix), art.opt-fbp.blur, imFormat)
end

MakePath('%spng',ParentPath);
prefix = sprintf('%spng/reco_full__',ParentPath);
postfix = '';
imFormat = 'png';
% Recos
WriteImage( sprintf('%sfbp_full%s',prefix,postfix), fbp.full, imFormat)
WriteImage( sprintf('%sfbp_red%s',prefix,postfix), fbp.red, imFormat)
WriteImage( sprintf('%sfbp_blur%s',prefix,postfix), fbp.blur, imFormat)
WriteImage( sprintf('%s%s_last__iter%04u%s',prefix,recType,NumIter,postfix), art.last, imFormat)
% Differences
WriteImage(sprintf('%sdiff__fbp_red-fbp_full%s',prefix,postfix), fbp.red-fbp.full, imFormat)
WriteImage( sprintf('%sdiff__fbp_blur-fbp_full%s',prefix,postfix), fbp.blur-fbp.full, imFormat)
WriteImage( sprintf('%sdiff__%s_last-fbp_full__iter%04u%s',prefix,recType,NumIter,postfix), art.last-fbp.full, imFormat)
WriteImage( sprintf('%sdiff__%s_last-fbp_blur__iter%04u%s',prefix,recType,NumIter,postfix), art.last-fbp.blur, imFormat)
if fbpErrMinPos < NumIter
    WriteImage( sprintf('%s%s_opt__iter%04u%s',prefix,recType,fbpErrMinPos,postfix), art.opt, imFormat)
    WriteImage( sprintf('%sdiff__%s_opt-fbp_full__iter%04u%s',prefix,recType,fbpErrMinPos,postfix), art.opt-fbp.full, imFormat)
    WriteImage( sprintf('%sdiff__%s_opt-fbp_blur__iter%04u%s',prefix,recType,fbpErrMinPos,postfix), art.opt-fbp.blur, imFormat)
end

%% ROI
% Crop recos
fbp.full = fbp.full(x,x);
fbp.red  = fbp.red(x,x);
fbp.blur = fbp.blur(x,x);
art.last  = art.last(x,x);
if fbpErrMinPos < NumIter
    art.opt  = art.opt(x,x);
end
% 
prefix = sprintf('%sreco_roi__',ParentPath);
postfix = '';
imFormat = 'png';
% Recos
WriteImage( sprintf('%sfbp_full%s',prefix,postfix), fbp.full, imFormat)
WriteImage( sprintf('%sfbp_red%s',prefix,postfix), fbp.red, imFormat)
WriteImage( sprintf('%sfbp_blur%s',prefix,postfix), fbp.blur, imFormat)
WriteImage( sprintf('%s%s_last__iter%04u%s',prefix,recType,NumIter,postfix), art.last, imFormat)
% Differences
WriteImage(sprintf('%sdiff__fbp_red-fbp_full%s',prefix,postfix), fbp.red-fbp.full, imFormat)
WriteImage( sprintf('%sdiff__fbp_blur-fbp_full%s',prefix,postfix), fbp.blur-fbp.full, imFormat)
WriteImage( sprintf('%sdiff__%s_last-fbp_full__iter%04u%s',prefix,recType,NumIter,postfix), art.last-fbp.full, imFormat)
WriteImage( sprintf('%sdiff__%s_last-fbp_blur__iter%04u%s',prefix,recType,NumIter,postfix), art.last-fbp.blur, imFormat)
if fbpErrMinPos < NumIter
    WriteImage( sprintf('%s%s_opt__iter%04u%s',prefix,recType,fbpErrMinPos,postfix), art.opt, imFormat)
    WriteImage( sprintf('%sdiff__%s_opt-fbp_full__iter%04u%s',prefix,recType,fbpErrMinPos,postfix), art.opt-fbp.full, imFormat)
    WriteImage( sprintf('%sdiff__%s_opt-fbp_blur__iter%04u%s',prefix,recType,fbpErrMinPos,postfix), art.opt-fbp.blur, imFormat)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if showFigs
    itool(fbp.full,'FBP full')
    itool(fbp.red,'FBP reduced')
    itool(fbp.blur,'FBP blurred')
    itool(art.last,'ART last iteration')
    if exist('art.opt','var')
        itool(art.opt,'ART optimal iteration (minimum of error art - fbp_full)')
    end
end