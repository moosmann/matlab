function cf_ = fit_piezo(test1,test2,amp,ampoff,phoff,period)
%FIT_PIEZO    Create plot of datasets and fits
%   FIT_PIEZO(TEST1,TEST2)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1

 
% Data from dataset "test2 vs. test1":
%    X = test1:
%    Y = test2:
%    Unweighted
%
% This function was automatically generated on 03-Mar-2011 14:36:06

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[449 168 621 429]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_); hold on;

 
% --- Plot data originally in dataset "test2 vs. test1"
test1 = test1(:);
test2 = test2(:);
h_ = line(test1,test2,'Parent',ax_,'Color',[0.333333 0 0.666667],...
     'LineStyle','none', 'LineWidth',1,...
     'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(test1));
xlim_(2) = max(xlim_(2),max(test1));
legh_(end+1) = h_;
legt_{end+1} = 'test2 vs. test1';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end


% --- Create fit "fit 1"
fo_ = fitoptions('method','NonlinearLeastSquares','Robust','On','Lower',[amp(2) period(2) phoff(2) ampoff(2)],'Upper',[amp(3) period(3) phoff(3) ampoff(3)]);
ok_ = isfinite(test1) & isfinite(test2);
st_ = [amp(1) period(1) phoff(1) ampoff(1) ];
set(fo_,'Startpoint',st_);
ft_ = fittype('a*sin(2*pi*(b*x+c))+d',...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'a', 'b', 'c', 'd'});

% Fit this model using new data
cf_ = fit(test1(ok_),test2(ok_),ft_,fo_);

% Or use coefficients from the original fit:
if 0
   cv_ = { 2, 0.06666666666667, 6.757183705948e-17, 3};
   cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
     'LineStyle','-', 'LineWidth',2,...
     'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 1';

% Done plotting data and fits.  Now finish up loose ends.
hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
xlabel(ax_,'');               % remove x label
ylabel(ax_,'');               % remove y label
