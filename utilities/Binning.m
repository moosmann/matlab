function array_binned = Binning( array, bin )
% Integer binning of 2D or 3D arrays. Array is cropped before binning
% such that mod(size(array), bin) = 0. Output array if of class 'single'.
%
% ARGUMENTS
% array : 2D or 3D array
% bin : integer scalar or 2- or 3- component vector. Default: 2. size of
%   bin for  each dimension. if scalar then bin is the same for each
%   dimension.
%
% OUTPUT
% array_binned : array, single.
%
% Written by Julian Moosmann.
%
% array_binned = Binning( array, bin)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    bin = 2;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if prod(bin) == 1 || isempty( bin ) || prod(bin) == 0
    array_binned = single( array );
    return
else
    switch ndims( array )
        
        case 2
            
            if isscalar( bin )
                bin1 = bin;
                bin2 = bin;
            else
                bin1 = bin(1);
                bin2 = bin(2);
            end
            
            [dim1,dim2] = size( array );
            
            % last relevant pixel
            last1 = dim1 - mod( dim1, bin1 );
            last2 = dim2 - mod( dim2, bin2 );
            
            %             array_binned = 0;
            %             array_binned = array_binned + single( array(1:2:last1,1:2:last2) );
            %             array_binned = array_binned + single( array(2:2:last1,1:2:last2) );
            %             array_binned = array_binned + single( array(1:2:last1,2:2:last2) );
            %             array_binned = array_binned + single( array(2:2:last1,2:2:last2) );
            %
            % Generic 2D bin functions
            %f = @(n,m)  cast( array(n:bin1:last1,m:bin2:last2), 'single' );
            %f = @(n,m) single( array(n:bin1:last1,m:bin2:last2) );
            f = @(n,m) array(n:bin1:last1,m:bin2:last2);
            
            % Binning
            array_binned = single( 0 );
            for xx = 1:bin1
                for yy = 1:bin2
                    %array_binned = array_binned + f(xx,yy);
                    array_binned = array_binned + single( f(xx,yy) );
                end
            end
            
        case 3
            
            if isscalar( bin )
                bin1 = bin;
                bin2 = bin;
                bin3 = bin;
            else
                bin1 = bin(1);
                bin2 = bin(2);
                bin3 = bin(3);
            end
            
            [dim1,dim2,dim3] = size( array );
            
            % crop to appropriate number of pixels
            last1 = dim1 - mod( dim1, bin1);
            last2 = dim2 - mod( dim2, bin2);
            last3 = dim3 - mod( dim3, bin3);
            
            % Generic 3D bin functions
            f = @(n,m,k) array(n:bin1:last1,m:bin2:last2,k:bin3:last3);
            
            array_binned = single( 0 );
            for xx = 1:bin1
                for yy = 1:bin2
                    for zz = 1:bin3
                        array_binned = array_binned + single( f(xx,yy,zz) );
                    end
                end
            end
            
        otherwise
            error( 'Binning of %u-D array not implemented.', ndims( array ) )
    end
    
end

