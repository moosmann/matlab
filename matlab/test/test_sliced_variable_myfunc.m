% function [a,s1,s2] = test_sliced_variable_myfunc(a, b, c)
% s1 = size(a,1);
% s2 = size(a,2);

% if b
%    fprintf( '\n Start loop' ) 
%     for n = 1:size(a,3)
%         
%         im = a(:,:,n);
%         im = im + 1;
%         im = padarray(im, size(im), 'post');
%         
%         a(:,:,n) = im(1:s1,1:s2);
%         if n ==1
%             fprintf( '\n First iteration finished' )
%         end
%         
%     end
% else
%     fprintf('\n c: %f' , c )
% end
% fprintf( '\n End of loop' )
% fprintf( '\n' )
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function proj = test_sliced_variable_myfunc( proj, phase_retrieval )

method = phase_retrieval.method;
padding = 1;
phase_filter = 1;
im_shape(1) = size(proj,1);
im_shape(2) = size(proj,2);


if strcmpi( method, 'tieNLO_Schwinger' )

     xi = 1;
     eta = 1;
    parfor (nn = 1:size(proj, 3), 4)
        
        g = SubtractMean( proj(:,:,nn) );
        g = gpuArray( g );
        
        g = padarray( g, padding * im_shape, 'symmetric', 'post' );
      
        %FT of g                        
        g_fts        = repmat( fft2( g ), [1, 1, sn] );
        
        % LO: zeroth order result        
        phi0 = 1 / ( 2 * pi ) * sum( ( ifft2( selap .* g_fts ) ), 3) / sn;
        phi0 = real(phi0);
        
        % NLO: first order result, first correction of three        
        phi11 = -1 / ( 4 * pi ) * sum( ifft2( selap .* repmat( fft2( g.^2 ),[1,1,sn] ) ), 3 ) / sn;
        % NLO: first order result, second correction
        phi12 = -1/(4*pi)*sum(real(ifft2(selap.*repmat(fft2( ...
            sum((ifft2( xi.*selap.*g_fts)),3)/sn.*ifft2( xi(:,:,1).*g_fts(:,:,1)) ...
            +sum((ifft2(eta.*selap.*g_fts)),3)/sn.*ifft2(eta(:,:,1).*g_fts(:,:,1)) ...
            ),[1,1,sn]))),3)/sn;
        phi13 = -1/(8*pi)*((sum((ifft2(xi.*selap.*g_fts)),3)/sn).^2 ...
            +(sum((ifft2(eta.*selap.*g_fts)),3)/sn).^2);        
        phi11 = real(phi11);
        phi13 = real(phi13);
        
        % LO + NLO
        g = phi0 + phi11 + phi12 + phi13;
        g = gather( g(1:im_shape(1),1:im_shape(2)) );
                
        proj(:,:,nn) = -g;
    end
else
    if phase_retrieval.use_parpool
        fprintf( 'using parpool' )
        parfor nn = 1:size(proj, 3)
            im = proj(:,:,nn);
            im = padarray( im, padding * im_shape, 'symmetric', 'post' );
            im = fft2( im );
            im = phase_filter .* im ;
            im = ifft2( im );
            im = -real( im );
            im = im(1:im_shape(1), 1:im_shape(2));
            proj(:,:,nn) = im;
        end
    else
        fprintf( 'without parpool' )
        for nn = 1:size(proj, 3)
            im = proj(:,:,nn);
            im = padarray( im, padding * im_shape, 'symmetric', 'post' );
            im = fft2( im );
            im = phase_filter .* im ;
            im = ifft2( im );
            im = -real( im );
            im = im(1:im_shape(1), 1:im_shape(2));
            proj(:,:,nn) = im;
        end
    end
end
