h_range = [12 35]; %mm
v_range = [40 55]; % mm

n = 1;
s(n).z = 50;
s(n).y = 20;
s(n).x = 10;

n = 2;
s(n).z = 55;
s(n).y = 15;
s(n).x = 10;

n = 3;
s(n).z = 40;
s(n).y = 12;
s(n).x = 10;

res = 5; %micron
exp_time = 100;%ms
fov = [7 2]; % mm
h = fov(1);
v = fov(2);

for n = 1:3
    x = s(n).x;
    y = s(n).y;
    z = s(n).z;

    num_h = ceil(max([x y]) / h);
    num_v = ceil(z / v);

    num_proj_nyquist = ceil(sqrt(2)*max([x y]) * 1000 * pi/2 / res / 100) * 100;
    scan_time_fov = num_proj_nyquist * exp_time / 1000;
    scan_time_sample_min = min(num_h) * min(num_v) * scan_time_fov / 60;

    fprintf('\n n: %u',n)
    fprintf('\n%20s: %u %u','step horizontal',num_h)
    fprintf('\n%20s: %u %u','step verical',num_v)
    fprintf('\n%20s: %u','# projectins Nyquist',num_proj_nyquist)
    fprintf('\n%20s: %.0f s = %.1f min','scan time fov',scan_time_fov,scan_time_fov/60)
    fprintf('\n%20s: %.1f min = %.1f h = %0.1f shifts = %0.1f d','scan time sample min',1./[1 60 60*8 60*24]*scan_time_sample_min)

end



fprintf('\n')

