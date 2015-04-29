y1 = afintr0p5;
y2 = afintr1p5;
y3 = afintr2p5;

x1 = x0p5;
x2 = x1p5;
x3 = x2p5;

ynorm = @(y,nn) y(:,nn) / max( y(:,nn) ) ;
ynormroi = @(y,nn,x) y(x,nn) / max( y(x,nn) ) ;


for nn=1:size(y1,2),
    plot( [ 1+ynormroi(y1,nn,x1) 0.5+ynormroi(y2,nn,x2) ynormroi(y3,nn,x3)] ,'.')
%     subplot(3,1,1)
%     plot( ynormroi(y1,nn,x1) )
%     subplot(3,1,2)
%     plot( ynormroi(y2,nn,x2) )
%     subplot(3,1,3)
%     plot( ynormroi(y3,nn,x3) )
%     axis tight
%     pause(0.3)
end