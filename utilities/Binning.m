function array = Binning(array, bin)
% Integer binning of 2D or 3D arrays. Array is cropped before binning
% such that mod(size(array), bin) = 0.
%
% array : 2D or 3D array
% bin: integer. Default: 2. size of bin
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

if isempty( bin ) || bin == 0 || bin == 1
    return
else
    switch ndims( array )
        case 2
            [x,y] = size( array );
            
            % last relevant pixel
            xl = x - mod( x, bin );
            yl = y - mod( y, bin );
            
            % Generic 2D bin functions
            f = @(n,m) array(n:bin:xl,m:bin:yl);
            
            % Binning
            tmp = 0;
            for mm = 1:bin
                for nn = 1:bin
                    tmp = tmp + f(mm,nn);
                end
            end
            array = tmp;
            
        case 3
            [x,y,z] = size( array );
            
            % crop to appropriate number of pixels
            xl = x - mod( x, bin);
            yl = y - mod( y, bin);
            zl = z - mod( z, bin);
            
            % Generic 3D bin functions
            f = @(n,m,k) array(n:bin:xl,m:bin:yl,k:bin:zl);
            
            % Binning
            switch bin
                case 2
                    tmp = f(1,1,1) + f(1,1,2) + f(1,2,1) + f(1,2,2) + f(2,1,1) + f(2,1,2) + f(2,2,1) + f(2,2,2);                
                otherwise
                    tmp = 0;
                    for ll = 1:bin
                        for mm = 1:bin
                            for nn = 1:bin
                                tmp = tmp + f(ll,mm,nn);
                            end
                        end
                    end                                        
            end
            array = tmp;
        
        otherwise
            error( 'Binning of %u-D array not implemented.', ndims( array ) )
    end
    
end

