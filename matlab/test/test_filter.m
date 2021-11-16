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
