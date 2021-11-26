p = '/asap3/petra3/gpfs/p05/2021/data/11010737/processed/02_AC11/reco_phase/tie_regPar1p80/float_rawBin2';
f = sprintf( '%s/reco_02_AC11_000300.tif', p );
im = imread( f );
imf0 = imsharpen( im );

uradius = 3;
uamount = 1;
uthresh = 0.0;
imf1 = imsharpen( im, 'Radius', uradius, 'Amount', uamount, 'Threshold', uthresh);

imd0 = -(imf0 - im);
imd1 = -(imf1 - im);

[s1, s2] = size(im);
s = 0.47;
s1r = ceil(s * s1):floor((1-s)*s1);
s2r = ceil(s * s2):floor((1-s)*s2);
s2rl = ceil(s * s2):ceil(0.5 * s2);
s2rr = floor(0.5 * s2):floor((1-s)*s2);

sp1 = 1;
sp2 = 2;

% subplot(sp1,sp2,1)
% imsc(im(s1r,s2r))
% xticklabels({});yticklabels({});axis off equal fill tight image
%
% subplot(sp1,sp2,2)
% imsc(imf(s1r,s2r))
% xticklabels({});yticklabels({});axis off equal fill tight image

%subplot(sp1,sp2,1)
%imsc(cat(2, im(s1r,s2rl), imf(s1r,s2rr) ))
mon = @(im,imf0,imd0) cat(2, normat( cat(2, im(s1r,s2r), imf0(s1r,s2r) ) ), normat( imd0(s1r,s2r) ) ) ;
imsc( cat(1, mon(im,imf0,imd0), mon(im,imf1,imd1) ) )
xticklabels({});yticklabels({});axis off equal fill tight image


% subplot(sp1,sp2,2)
% imsc(imd(s1r,s2r))
% xticklabels({});yticklabels({});axis off equal fill tight image

% subplot(sp1,sp2,4:6)
% imsc(imd(s1r,s2r))
% xticklabels({});yticklabels({});axis off equal fill tight image

%% EMBL background filter
%p = '/asap3/petra3/gpfs/p05/2020/data/11009667/raw/embl_060_200924_dist_1400_zshift_0p2_010';

fn = '/asap3/petra3/gpfs/p05/2020/data/11009667/processed/images/flat_fields_dark_corrected.mat';
load( fn )
%%

im = rot90(flat(:,:,3));

rect1 = [1608         113         180         163];
rect2 = [1412 17 614 635];

[patch1] = imcrop(im, rect1);
[patch2] = im;%imcrop(im, rect2);

patchVar1 = std2(patch1)^2;
patchVar2 = std2(patch2)^2;

sf = 1;
DoS1 = sf*2*patchVar1;
DoS2 = sf*2*patchVar2;

imf1 = imbilatfilt(im,DoS1);
imf2 = imbilatfilt(im,DoS2);

if exist( 'h1' , 'var' ) && isvalid( h1 )
    figure(h1)
else
    h1 = figure( 'Name', 'bilateral filter' );
end

subplot(3,1,1)
imsc( im )
axis equal tight
drawrectangle( 'Position', rect )

subplot(3,1,2)
imsc( imf1 )
axis equal tight

subplot(3,1,3)
imsc( imf2 )
axis equal tight

if exist( 'h2' , 'var' ) && isvalid( h2 )
    figure(h2)
else
    h2 = figure( 'Name', 'bilateral filter difference' );
end

subplot(3,1,1)
imsc( imf1 - imf2 )
axis equal tight

subplot(3,1,2)
imsc( imf1 - im)
axis equal tight

subplot(3,1,3)
imsc( imf2 - im)
axis equal tight