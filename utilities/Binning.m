function array = Binning(array, bin)
% Integer binning of 2D or 3D arrays. Array is cropped before binning
% such that mod(size(array), bin) = 0.
%
% array : 2D or 3D array
% bin: integer scalar or 2- or 3- component vector. Default: 2. size of bin for
% each dimension. if scalar then bin is the same for each dimension
%
% Written by Julian Moosmann.
% Last modification 2017-10-09
%
% array = Binning(array, bin)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    bin = 2;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty( bin ) || prod(bin) == 0 || prod(bin) == 1
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
            
            % Generic 2D bin functions
            f = @(n,m) array(n:bin1:last1,m:bin2:last2);
            
            % Binning
            tmp = 0;
            for mm = 1:bin1
                for nn = 1:bin2
                    tmp = tmp + f(mm,nn);
                end
            end
            array = tmp;
            
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
            
            % Binning
%             switch bin1 
%                 case 2
%                     tmp = f(1,1,1) + f(1,1,2) + f(1,2,1) + f(1,2,2) + f(2,1,1) + f(2,1,2) + f(2,2,1) + f(2,2,2);                
%                 otherwise
                    tmp = 0;
                    for ll = 1:bin1
                        for mm = 1:bin2
                            for nn = 1:bin3
                                tmp = tmp + f(ll,mm,nn);
                            end
                        end
                    end                                        
%            end
            array = tmp;
        
        otherwise
            error( 'Binning of %u-D array not implemented.', ndims( array ) )
    end
    
end

