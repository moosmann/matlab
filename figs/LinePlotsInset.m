function LinePlotsInset(exact_phase,phi,distance,line_cut_at,plot_range_zoomed,plot_range)
    
    if nargin<4,line_cut_at=floor(size(exact_phase,2)/3);end;
    if nargin<6,plot_range=0;end;
    if nargin<5,plot_range_zoomed=0;end;

resolution = size(exact_phase,1);
prefix = '/home/moosmann/data/SiemensStar/paper_figures';
frame_title = sprintf(['Line cut at %g: Black=ExactPhase, Blue=Bronnikov, Red=Bronnikov+Correction'],line_cut_at);
if ~plot_range,
    dimx = size(exact_phase,1);
    xini = mod(dimx,100);
    x    = (1:dimx-xini);
else, %x = floor(size(exact_phase)/4)+plot_range;end;
      x = plot_range;end;
if line_cut_at==0,
    y = floor(size(exact_phase,2)/3);
else,
    y = line_cut_at;
end;
if plot_range_zoomed==0,
    xz = 234:244;
else,
    xz = plot_range_zoomed;
end;

%set(0,'defaulttextinterpreter','none');
f = figure('Name',frame_title);
set(f,'color','white');
axes('Box','on','XTick',[],'YTick',[]);
a1=axes('position',[.13 .11 .774 .8125]);
p1=plot(x,exact_phase(x,y),'black', ... 
     x,phi(x,y,1),'blue', ... 
     x,phi(x,y,1)+phi(x,y,2),'red');
%ylabel(sprintf('phase shift $\\phi(x,y=%g)-\\langle\\phi\\rangle$',line_cut_at));
%xlabel('pixel');
set(a1,'TickDir','out','Box','off');
%set(gca,'Box','off'),
axes('position',[.455 0.17 0.405 0.405],'Box','on','XTick',[],'YTick',[]);
a2in=axes('position',[.46 0.17 0.4 0.4]);
p2in=plot(xz,exact_phase(xz,y),'black', ... 
     xz,phi(xz,y,1),'blue', ... 
     xz,phi(xz,y,1)+phi(xz,y,2),'red');
%set(p2in,'Marker','.','MarkerSize',5.0);
set(a2in,'Box','off','YAxisLocation','right','TickDir','out');
%file_name = sprintf('%s/intensity_and_phase_maps',prefix);
file_name = sprintf('%s/LineCut_z%03gcm_res%04g',prefix,100*distance,resolution);
saveas(f,file_name,'epsc');
%laprint(f,file_name);





