tic
%parfor (nn =  1:100)
for nn =  1:100
    im = proj0(:,:,nn);
    im = NegLog(im, take_neg_log);
    im = padarray( im, padding * [proj_shape1 0 0], padding_method, 'post' );
    im = fft( im, [], 1);
    im = bsxfun(@times, im, filt);
    im = real( ifft( im, [], 1, 'symmetric') );
    im = im(1:proj_shape1,:,:);
    proj(:,:,nn) = im;
end
toc