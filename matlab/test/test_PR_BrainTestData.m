

p = '/beegfs/desy/user/moosmanj/BrainTestData';
if ~exist('s','var')
    s = read_images_to_stack(p);
end

[s1, s2, s3] = size(s);
method = 'tie';
edp = [1.2398e-06 2/pi 1];
reg_par = -20;
[phase_filter, str] = PhaseFilter(method,2*[s1 s2],edp,reg_par);
padding = 1;

sp = zeros(size(s),class(s));
for n=1:s3
    im = s(:,:,n);
    %im = SubtractMean(im);
    im = padarray( im, padding * [s1 s2], 'symmetric', 'post' );
    im = fft2( im );
    im = phase_filter .* im ;
    im = ifft2( im );
    im = real( im );
    sp(:,:,n) = im(1:s1, 1:s2);
    savepath = sprintf('%s/%s',p,str);
    CheckAndMakePath(savepath)
    fn = sprintf('%s/phase_%04u.tif',savepath,n);
    write32bitTIFfromSingle(fn,sp(:,:,n));
end