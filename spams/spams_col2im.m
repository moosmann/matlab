function [im, mask] = spams_col2im(mat,blockSize,imSize)
% Rearrange matrix which is made from image blocks of 'blockSize' back to
% image format.

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    blockSize = [8 8];
end
if nargin < 3
    imSize = [512 512];
end

%% Main 
m = blockSize(1);
n = blockSize(2);
M = imSize(1);
N = imSize(2);

im = zeros( M, N);
mask = zeros( M, N);
patch1 = ones(m,n);

% Add patches up
for ii = 1:M-m+1
    for jj = 1:N-n+1
        %im( ii+(0:m-1), jj+(0:n-1) ) = reshape(mat(:,ii + (M-n+1)*(jj-1)), [m n] ) ;
        patch = reshape(mat(:,ii + (M-n+1)*(jj-1)), [m n] ) ;
        im( ii+(0:m-1), jj+(0:n-1) ) = im( ii+(0:m-1), jj+(0:n-1) ) + patch;
        mask( ii+(0:m-1), jj+(0:n-1) ) = mask( ii+(0:m-1), jj+(0:n-1) ) + patch1;
    end
end

% Normalize image
im = im./mask;
