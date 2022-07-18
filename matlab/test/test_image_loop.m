x = round( [0.2 0.8] * (size(vol,1) - 1) + 1);
y = round( [0.2 0.8] * (size(vol,2) - 1) + 1);
xr = x(1):x(2);
yr = y(1):y(2);
tic;
fprintf(' \n  0 ')
im = vol(xr,yr,506+20*nn);
imi = uint8(2^8 * normat( im ));
h = imagesc(imi);
axis image tight
fprintf(' %10.2f',toc);
hold on;
for nn = 1:3
    fprintf('\n %2u',nn);
    %hold on;
    %imagesc('CData',vol(xr,yr,506+50*nn));
    im = vol(xr,yr,506+20*nn);
    imi = uint8(2^8 * normat( im ));
    %imagesc(imi)
    imagesc('CData', imi)
    %h.CData = vol(x,y,506+20*nn);
    %refreshdata
    drawnow;
    fprintf(' %10.2f',toc);
end
fprintf('\n')