function N = functest(N)

if nargin < 1
    N = 5;
end

disp(N);
% fprintf('\nInput argument: %g\n',N)
% N = 2*N;
% fprintf('\nDouble input argument: %g\n',N)

%% stupid
center = floor((N + 1)/2);
xleft = -center ;
x = (1:N) + xleft;
x = repmat(x, N, 1);

ytop = center - 1;
y = (N:-1:1).' - N + ytop;
y = repmat(y, 1, N);

disp(x),disp(y)

%% smart
x0=-floor((N-1)/2):floor(N/2);
[xx,yy]=meshgrid(x0,-x0);
disp(xx),disp(yy)

disp([isequal(x,xx),isequal(y,yy)])

%% check interpolation
theta=0:18:180;costheta = cos(theta);sintheta = sin(theta);i=1;t=x.*costheta(i)+y.*sintheta(i);t=t(:);t2=x(:).*costheta(i)+y(:).*sintheta(i);isequal(t,t2)