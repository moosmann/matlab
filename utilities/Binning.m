function array = Binning(array, bin)
% 2 x 2 or 4 x 4 binning of 2-D image. Images is cropped before binning
% such that mod(size(array), bin) = 0.
%
% array: 2D image to bin
% bin: bin size, 2 or 4. Default: 2
%
% Written by Julian Moosmann.
% Last modification 2016-12-06
%
% array = Binning(array, bin)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    bin = 2;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bin == 0 || bin == 1
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
            
            switch bin
                case 2
                    array = f(1,1) + f(1,2) + f(2,1) + f(2,2);
                case 3
                    array = f(1,1) + f(1,2) + f(1,3) + f(2,1) + f(2,2) + f(2,3) + f(3,1) + f(3,2) + f(3,3);                          
                case 4
                    array = f(1,1) + f(1,2) + f(1,3) + f(1,4) + f(2,1) + f(2,2) + f(2,3) + f(2,4) + f(3,1) + f(3,2) + f(3,3) + f(3,4) + f(4,1) + f(4,2) + f(4,3) + f(4,4);
                otherwise
                    error('Bin size %g not implemented.', bin)
            end
            
        case 3            
            [x,y,z] = size( array );
            
            % crop to appropriate number of pixels            
            xl = x - mod( x, bin);
            yl = y - mod( y, bin);
            zl = z - mod( z, bin);
                                    
            % Generic 3D bin functions            
            f = @(n,m,k) array(n:bin:xl,m:bin:yl,k:bin:zl);
            
            switch bin
                case 2
                    array = f(1,1,1) + f(1,1,2) + f(1,2,1) + f(1,2,2) + f(2,1,1) + f(2,1,2) + f(2,2,1) + f(2,2,2);
                case 3                    
                    error('Bin size %g not implemented.', bin)
                case 4
                    error('Bin size %g not implemented.', bin)                    
                otherwise
                    error('Bin size %g not implemented.', bin)
            end
            
        otherwise
            error( 'Binning of %u-D array not implemented.', ndims( array ) )
    end
    
end

