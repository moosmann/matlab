function [cbar,figure1,axes1] = savecolorbar(climits,filename)

% Create grayscale colorbar with limits [climits(1) climits(2)], or if
% input climits is an image climits = [min(climits(:) max(climits(:))] is
% used and save figure in eps format.
%
% Written by J. Moosmann, first version: 2014-06-24

if nargin < 1
    climits = [0.92 1.08];
end
if nargin < 2
    filename = 'colorbar';
end

%% Main
fontsize = 28;
if size(climits,1) > 1
    climits = [min(climits(:)) max(climits(:))];
end
    %,'FontSize',24);
% Create figure
figure1 = figure;
axes1 = axes('Visible','off','Parent',figure1,'CLim',climits,'FontSize',fontsize);
cbar = colorbar('peer',axes1,[0.1 0.13 0.3 0.77],'FontSize',fontsize);
colormap('gray');

set(figure1,'PaperPositionMode','Manual')
set(figure1,'PaperUnits','normalized','PaperPosition',[0 0 .3 1]);
%set(figure1,'Units','normalized','OuterPosition',[0.5 0.5 0.1 .5]);


set(axes1,'FontSize',fontsize);

% Save figure
saveas(figure1,[filename '.eps'])
%saveas(figure1,[filename '.png'])


