% Plot filter kernel.

dimx = 16;
dimy = 32;
[f1 f2]=freqspace([dimx dimy],'meshgrid');
Hd = zeros(dimx,dimy);
Hd(floor(dimx/4):end-ceil(dimx/4),floor(dimy/4):end-ceil(dimy/4))=1;

win = fspecial('gaussian',[dimx dimy],2);
win = win ./ max(win(:));

h = fwind2(Hd,win);

if 1,
    figure('Name','window'),colormap(jet(64)),mesh(f1,f2,Hd)
    figure('Name','Gaussian filter'),mesh(win)
    figure('Name','frequency response'),freqz2(h)
end
