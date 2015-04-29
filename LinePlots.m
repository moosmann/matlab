function LinePlots(exact_phase,phi,line_cut_at,plot_range,frame_title_prefix)
    
% LinePlots(exact_phase,phi,line_cut_at,plot_range,frame_title_prefix)
    
    if nargin<3,line_cut_at=floor(size(exact_phase,2)/3);end;
    if nargin<4,plot_range=0;end;
    if nargin<5,frame_title_prefix='Retrieved phase: ';end;
%frame_title = sprintf('%sLine cut at %g: Black=ExactPhase, Blue=Bronnikov, Red=Bronnikov+Correction',frame_title_prefix,line_cut_at);
frame_title = '';
if ~plot_range, x = (1:size(exact_phase,2));
else %x = floor(size(exact_phase)/4)+plot_range;end;
      x = plot_range;
end;
y = line_cut_at;

figure('color','white','Name',frame_title);
plot(x,exact_phase(x,y),'black', ... 
     x,phi(x,y,1),'blue', ...
    x,phi(x,y,1)+phi(x,y,2),'red'), ... 
    title(frame_title);
set(gca,'FontSize',16);


