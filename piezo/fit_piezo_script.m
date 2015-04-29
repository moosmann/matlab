%fit_piezo_script
gridy_delta = 1.448;
y1=squeeze(im(1024,1024,1:16));
y2=squeeze(im(1024+1,1024,1:16));
x=gridy_delta*(0:numel(y1)-1)';
plot(x,y1,x,y2);
cf1 = fit_piezo(x,y1,[1000 700 1800],[2800 2500 3800],[0 -16 16],1/gridy_delta*[1/8 1/16 1/4]);
cf2 = fit_piezo(x,y2,[1000 700 1800],[2800 2500 3800],[0 -16 16],1/gridy_delta*[1/8 1/16 1/4]);
cf = [cf1.b cf2.b]';
disp(cf);
gridy_cal = mean(cf);
disp(gridy_cal);
