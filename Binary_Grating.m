function Binary_Grating(period,displacement,m,n,sampling)
    
    k=2*pi*(1:sampling)/sampling/displacement;
    fg = (m*sin(k*period/2)./k).^2.*(1+(n/m)^2+n/m*sin(k*(m*period+ ...
                                                      displacement)));
    figure(1),plot(1:sampling,fg);
    