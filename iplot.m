function iplot(data,alongDim,pause_s)
% Play plot
%
% Written by Julian Moosmann, 2015-04-14

%% Default arguments
if nargin < 2
    alongDim = 2;
end
if nargin < 3
    pause_s = 0.1;
end

%% MAIN
%h = figure('Name','test');
set( gcf, 'Name', sprintf( '%s. Play plot', inputname(1)) );
colormap(gray)

switch alongDim
    case 1
        for nn = 1:size( data, 1)
            plot( data(nn,:) )
            axis tight;
            pause( pause_s );
        end
    case 2
        for nn = 1:size( data, 2)
            plot( data(:,nn) )
            axis tight;
            pause( pause_s );
        end
end
