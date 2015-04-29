function im = FilterArtefact(im,roi_x1_y1_x2_y2)
% Replace corrupted region by median filtered region

x1 = roi_x1_y1_x2_y2(1);
y1 = roi_x1_y1_x2_y2(2);
x2 = roi_x1_y1_x2_y2(3);
y2 = roi_x1_y1_x2_y2(4);
imroi = im(x1:x2,y1:y2);
imroi = medfilt2(imroi,[ceil((x2-x1)/2+1) ceil((y2-y1)/2+1)],'symmetric');
im(x1:x2,y1:y2) = imroi;