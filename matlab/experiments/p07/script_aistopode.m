p = '/asap3/petra3/gpfs/p07/2024/data/11020243/raw/003_aistopode__a/';

if 1
ref = 10;
d = dir([p 'tiff0000/*ref.tif']);
for n=1:numel(d)
    fn = [d(n).folder filesep d(n).name];
    im = single(imread(fn));
    im = FilterPixel(im);
    ref = ref + im;
end
ref = Binning(ref,2);
end

im1 = single(imread([p 'tiff0001/003_aistopode__a_005419_img.tif']));
im2 = single(imread([p 'tiff0003/003_aistopode__a_015418_img.tif']));

im1 = Binning(FilterPixel(im1),2);
im2 = Binning(FilterPixel(im2),2);



im1 = im1./ref;
im2 = im2./ref;

im1 = im1(200:end-200,1:500);
im2 = fliplr(im2(200:end-200,1:500));

%itool([im1 im2])

c = ImageCorrelation(im1,im2);
x = int16(c.shift1);
y = int16(c.shift2);

im1c = im1(1+x:end,1+y:end);
im2c = im2(1:end-x,1:end-y);

cc = ImageCorrelation(im1c,im2c);

nimplay(cat(3,im1c',im2c'))