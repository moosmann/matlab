function plotbrontheta(stack,alpha,iterations,dp,print);
%error measure over all slice (angles)
inpnam=inputname(1);
for ii=1:dp:360;
  x(ii)=ii/2;
  [mean0(ii),mean(ii),diff0(ii),diff(ii)]=recbron(stack,ii,1,alpha,iterations,0);
end;

int=1:dp:360;
%figure('Name','Phantom: blue: mean(|phi_bron|), read: mean(|phi_better|)'), ... 
 %   plot(x(int),mean0(int),'blue',x(int),mean(int),'red');
figure('Name',sprintf('Phantom: PSI: dist%s blue-bro, red-brocor.eps',inpnam(5:6))), ... 
    plot(x(int),diff0(int),'blue',x(int),diff(int),'red'),
axis([0 180 2*10^-3 5.5*10^-3]);
if print
saveas(gcf,sprintf('phantom_PSI_d%s_blue_bro_red_brocor.eps',inpnam(5:6)),'psc2');
end;