function LinePlots2(phase,linTIE,CTF,x,y,FontSize)
    
    if nargin<4
        x=532;
    end
    if nargin<5
        y=450:600;
    end
    if nargin<6
        FontSize = 32;
    end

figure;
plot(y,linTIE(y,x,1),'blue', ... 
     y,CTF(y,x,1),'green', ...
     y,phase(y,x),'red');
 set(gca,'FontSize',FontSize);
