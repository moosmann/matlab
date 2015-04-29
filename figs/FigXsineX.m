function FigXsineX
% Plot x, sin(x), (x-sin(x))/x, and (x-sin(x))/sin(x).

x=(0:1/12:.75*pi);
figure,plot(x,sin(x),'black',x,x,'blue',x,(x-sin(x))./x,'yellow',x,(x-sin(x))./sin(x),'magenta');
set(gca,'XTick',0:pi/4:pi)
set(gca,'XTickLabel',{'0','pi/4','pi','3pi/4','pi'})
for x=[pi/12 pi/6 pi/4 pi/2]
       text(x,(x-sin(x))/x,['(' num2str(x,'%.2g') ',' num2str((x-sin(x))/x,'%3.2g') ...
                      ',' num2str((x-sin(x))/sin(x),'%3.2g') ')']);
end;
legend('x','sin(x)', '(x-sin(x))/x','(x-sin(x))/sin(x)','Location','NorthWest');