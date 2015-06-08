function iplot(data,alongDim,pause_s)
% Play plot
%
% Written by Julian Moosmann, 2015-04-14

%% Default arguments
if nargin < 2
    alongDim = 2;
end
if nargin < 3
    pause_s = 0;
end

%% MAIN
%h = figure('Name','test');
set( gcf, 'Name', sprintf( '%s. Play plot', inputname(1)) );
colormap(gray)



switch alongDim
    case 1
        nmax = size( data, 1);
        if pause_s == 0
            pause_s = 10 / nmax;
        end
        for nn = 1:nmax
            plot( data(nn,:), '.')
            legend(sprintf('%u of %u', nn, nmax))
            axis tight;
            pause( pause_s );
        end
    case 2
        nmax = size( data, 2);
        if pause_s == 0
            pause_s = 10 / nmax;
        end
        for nn = 1:nmax
            plot( data(:,nn), '.' )
            legend(sprintf('%u of %u', nn, nmax))
            axis tight;
            pause( pause_s );
        end
end
