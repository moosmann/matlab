p = '/asap3/petra3/gpfs/p05/2020/data/11009117/raw/images/';
darks = read_images_to_stack(p,1,'dark*.tif');
dark = median( dark, 3 );
domain( dark );

im1 = imread([p 'raw_image_30keV_07082020_after-fullfield_dosimetry.tif']) - dark;
im2 = imread([p 'raw_image_MSC_aligned_30keV.tif']) - dark;

domain( im1)
domain(im2)

im1 = FilterPixel(im1);
im2 = FilterPixel(im2);


domain( im1)
domain(im2)