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





if 0
    clear all
    
    v = round( 100 * phantom3d( 'modified shepp-logan', 512 ) );
    %v = zeros( [10 10 10] );
    %v( 2:8, 2:8, 2:8 ) = 20;
    %v( 4:6, 4:6, 4:6 ) = 30;
    m = zeros( size( v ) );
    
    [sx, sy, sz] = size( v );
    
    s = v( :,:, round( sz / 2 ) );
    
    
    %itool( s )
    l1 = 20;
    l2 = 30;
    
    c = 0;
    tic
    for x = 1:sx
        for y = 1:sy
            for z = 1:sz
                if isequal( v(x,y,z), l1 )
                    if x+1 < sx && isequal( v(x+1,y,z), l2 )
                        c = c + 1;
                        m(x,y,z) = m(x,y,z)+1;
                    end
                    if y+1 < sy && isequal( v(x,y+1,z), l2 )
                        c = c + 1;
                        m(x,y,z) = m(x,y,z)+1;
                    end
                    if z+1 < sz && isequal( v(x,y,z+1), l2 )
                        c = c + 1;
                        m(x,y,z) = m(x,y,z)+1;
                    end
                end
                
                if isequal( v(x,y,z), l2 )
                    if x+1 < sx && isequal( v(x+1,y,z), l1 )
                        c = c + 1;
                        m(x+1,y,z) = m(x+1,y,z)+1;
                    end
                    if y+1 < sy && isequal( v(x,y+1,z), l1 )
                        c = c + 1;
                        m(x,y+1,z) = m(x,y+1,z)+1;
                    end
                    if z+1 < sz && isequal( v(x,y,z+1), l1 )
                        c = c + 1;
                        m(x,y,z+1) = m(x,y,z+1)+1;
                    end
                end
                
                
            end
        end
    end
    t = toc;
    
    fprintf( '\n contact voxels : %u', c )
    fprintf( '\n time elapsed : %f', t )
    nimplay( v )
    nimplay( m )
    
end