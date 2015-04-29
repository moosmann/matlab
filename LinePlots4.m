function LinePlots4(exact_phase,phi1,phi2,line_cut_at,plot_range)
    
% LinePlots(exact_phase,phi,line_cut_at,plot_range,frame_title_prefix)
    
    if nargin<4,line_cut_at=floor(size(exact_phase,2)/3);end;
    if nargin<5,plot_range=0;end;

if ~plot_range, x = (1:size(exact_phase,2));
else %x = floor(size(exact_phase)/4)+plot_range;end;
      x = plot_range;
end;
y = line_cut_at;

figure('color','white');
plot(x,exact_phase(x,y),'black', ... 
     x,phi1(x,y),'blue', ...
    x,phi2(x,y),'red'),
set(gca,'FontSize',16);


