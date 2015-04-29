ca

n = 2048;

rl = (iradonDesignFilter('ram-lak',n,1))';

x = [ (-n+1):1:-1  0:1:n];


rl0 =fftshift( 1 / n * abs(x) );
disp(rl([1, n-1, n+1, n]))
disp(rl0([1, n-1, n+1, n]))
x = fftshift(x);
%plot(x,[rl0 ; rl]')
drl = rl-rl0;
xx = 0+(1:10);
plot(x(xx),drl(xx))
plot(x(xx),[rl(xx)-rl0(xx)],'.')