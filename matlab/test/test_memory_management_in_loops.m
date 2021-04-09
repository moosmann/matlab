% Test influence of memory allocation and variable management in loops.

nump = 30;
num_it = 109;

% fprintf( ' 2' )
% tic
% parfor (nn =  1:num_it, nump)
% %for nn =  1:num_it
%     im = proj0(:,:,nn); 
%     im = NegLog(im, take_neg_log);
%     imp = padarray( im, padding * [proj_shape1 0 0], padding_method, 'post' );
%     impc = fft( imp, [], 1);
%     impc = bsxfun(@times, impc, filt);
%     imp = ifft( impc, [], 1, 'symmetric');
%     imp = real( imp );
%     im = imp(1:proj_shape1,:,:);
%     proj(:,:,nn) = im;
% end
% t2 = toc;

tic
parfor (nn =  1:num_it, nump)
%for nn =  1:num_it
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

% 
% fprintf( ' 1' )
% tic
% parfor (nn =  1:num_it, nump)
% %for nn =  1:num_it
%     im = proj0(:,:,nn); 
%     tmp = NegLog(im, take_neg_log);
%     im = padarray( tmp, padding * [proj_shape1 0 0], padding_method, 'post' );
%     tmp = fft( im, [], 1);
%     im = bsxfun(@times, tmp, filt);
%     tmp = real( ifft( im, [], 1, 'symmetric') );
%     im = tmp(1:proj_shape1,:,:);
%     proj(:,:,nn) = im;
% end
% t1 = toc;



