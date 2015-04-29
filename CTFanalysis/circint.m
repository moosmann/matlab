%% Parameters.
tt = zeros( numel(theta), size(afint,3));
for jj = 1:size(afint,3)
thetaoffset = 1*0.05;
thetastep = 1/2/250;
theta = 0+thetaoffset:thetastep:pi/2-thetaoffset;
theta = theta(:);
theta = cat(1,theta,theta+pi/2,theta+pi,theta+3*pi/2);

%% Programme.

for rr = 511:-1:1
    for ii = length(theta):-1:1
        tt(ii,jj) = afint(round(512+rr*cos(theta(ii))),round(512+rr*sin(theta(ii))),jj);
    end
    afintLine(rr,jj) = sum(tt(:,jj))/length(theta);
end

%whos afintLine
x = 1:511;
%figure
plot(x,squeeze(afintLine(:,jj))),axis tight,pause(0.1)
end