function write_multiTiff_sequence( folder)

if nargin < 1
    folder = '/asap3/petra3/gpfs/p05/2021/data/11008741/processed/syn0168_102L_Mg5Gd_2d_0*';
end
if nargin < 2
    subfolder = 'dual_denoise';
end
if nargin < 3
    outpath = '/asap3/petra3/gpfs/p05/2021/data/11008741/processed/syn0168_102L_Mg5Gd_2d';
end
if nargin < 4
    outname = 'vol_';
end

CheckAndMakePath( outpath )

d = dir( folder );

for n = 11:numel(d)
    p = [d(n).folder filesep d(n).name filesep subfolder];
    v = read_images_to_stack( p, 1, '*.tif', [], 1, 1);
    
    fn = sprintf( '%s/%s%04u.am', outpath, outname, n);
    write_HxUniformScalarField(fn, v, 'single' )
    
%     fn = sprintf( '%s/%s%04u.tif', outpath, outname, n);
%     for s = 1:size(v, 3)
%         if s ~= 1            
%             write32bitTIFfromSingle( fn, v(:,:,s), 'a')
%         else
%             write32bitTIFfromSingle( fn, v(:,:,s), 'w')
%         end
%     end

end