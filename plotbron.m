function plotbron(proj,slice,alphamax,iterations,pts,cut1,cut2,fig);
%print figures of Bronnikov phantom
%print:sum_i |phi0_i|(alpha),sum_i |phi_i|(alpha),sum_i (|phi0_i|-|phi_i|)(alpha)

for ii=0:pts;
  x(ii+1)=ii/pts*alphamax;
  [mean0(ii+1),mean(ii+1),diff0(ii+1),diff(ii+1)] = recbron(proj,slice,1,x(ii+1),iterations,0);
end;

%figure: mean0,mean1
if fig==1
int=1:pts+1;
figure('Name','Phantom: blue:mean_bron, red:mean_we'), ... 
    plot(x(int),mean0(int),'blue',x(int),mean(int),'red');
cd ~/data/phantom/pics;
inpnam=inputname(1);
saveas(gcf,sprintf('phantom_PHI_d%s_th%03u_blue_bro_red_brocor.eps',inpnam(5:6),slice/2),'psc2');
end;
%figure: mean0-mean1
%int=1:pts+1;
%figure('Name','Phantom: mean_bron - mean_we'), ... 
%    plot(x(int),mean0(int)-mean(int),'blue');
%figure: zooming in mean0-mean1
%if cut1~=0 & cut2~=0
%int=1+(floor(cut1*pts/alphamax):floor(cut2*pts/alphamax));
%figure('Name','Phantom: mean_bron - mean_we, zoomed'), ... 
%    plot(x(int),mean0(int)-mean(int),'blue');
%end;


%figure: diff(phi_exact-phi0),diff(phi_exact-phi)
if (fig==2 || fig==3)
int=1:pts+1;
figure('Name','Phantom: blue: mean(|phi_exact - phi_bron|), red:mean(|phi_exact - phi_we|)'), ... 
       plot(x(int),diff0(int),'blue',x(int),diff(int),'red');
end;
%figure: zooming in of figure above
if cut1~=0 & cut2~=0
int=1+(floor(cut1*pts/alphamax):floor(cut2*pts/alphamax));
figure('Name','Phantom: blue:mean(|phi_exact - phi_bron|), red:mean(|phi_exact - phi_we|), zoomed'), ...
    plot(x(int),diff0(int),'blue',x(int),diff(int),'red');
end;




