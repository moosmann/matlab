% Compute blurring due to finite source size at detector position.

% Sample-detector distance in mm
z = (75+(0:1000/3:1000));
% Source size in mm
ds = 0.25;
% Source sample distance in mm: distance to Be window + Be-sample stage
% disance
s = (27500+1000);
% Blurring in mm
b = ds/s*z;
fprintf('\nSample-detector distance: %s mm\n',mat2str(z,4))
fprintf('Sourece size: %g mm\n',ds)
fprintf('Source-sample distance: %g m\n',s/1000)
fprintf('Blurring at detector position: %s mu\n\n',mat2str(b*1000,2))