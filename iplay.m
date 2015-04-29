function iplay(vol,along_dim,transpose_displayed_image,pause_s)
% Play volume 'vol' along dimension 'along_dim' using imagesc.
%
% Written by Julian Moosmann, 2015-03-24

%% Default arguments
if nargin < 2
    along_dim = 3;
end
if nargin < 3
    transpose_displayed_image = 0;
end
if nargin < 4
    pause_s = 0.05;
end

%% MAIN
%h = figure('Name','test');
set( gcf, 'Name', sprintf( '%s, %ux%ux%u, playing along dim %u', inputname(1), size(vol), along_dim ) );
colormap(gray)

switch along_dim
    case 1
        for nn = 1:size( vol, along_dim)
            if ~transpose_displayed_image
                imagesc( squeeze( vol(nn,:,:) ) );
            else
                imagesc( squeeze( vol(nn,:,:) )' );
            end
            axis equal tight off;
            pause( pause_s );
        end
    case 2
        for nn = 1:size( vol, along_dim)
            if ~transpose_displayed_image
                imagesc( squeeze( vol(:,nn,:) ) );
            else
                imagesc( squeeze( vol(:,nn,:) )' );
            end
            axis equal tight off;
            pause( pause_s );
        end
    case 3
        for nn = 1:size( vol, along_dim)
            if ~transpose_displayed_image                
                imagesc( squeeze( vol(:,:,nn) ) );
            else
                imagesc( squeeze( vol(:,:,nn) )' );
            end
            axis equal tight off;
            pause( pause_s );
        end
end

