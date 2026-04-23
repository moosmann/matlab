E_keV = 18;
l = E_to_lambda(E_keV*1000);
p = 0.65e-6;
z = 44e-3;
z = 304e-3;
dob = 0.99/0.11

F = p^2/l/z;

q = FrequencyVector(1024);
arg = pi*l*z/p^2/4 *abs(q).^2;
H = cos(arg) + dob*sin(arg);
plot(H)

