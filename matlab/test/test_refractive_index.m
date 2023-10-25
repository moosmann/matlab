clear n
n.h2o(1).energy_ev = 20000;
n.h2o(1).delta = 5.76620721E-07;
n.h2o(1).beta = 3.46201373E-10;
n.h2o(1).density_gcm3 = 1;

n.h2o(2).energy_ev = 30000;
n.h2o(2).delta = 2.56114134E-07;
n.h2o(2).beta = 1.06431752E-10;

n.c(3).density_gcm3 = 2.2;
n.c(3).energy_ev = 20000;
n.c(3).delta = 1.14115369E-06;
n.c(3).beta = 3.96117222E-10;

n.c(4).energy_ev = 30000;
n.c(4).delta = 5.0698884E-07;
n.c(4).beta = 1.57106508E-10;
    

z = 0.1e-3; %m
nn = n.c(3);
energy = nn.energy_ev;
lambda = ev_to_lambda(energy);
delta = nn.delta;
beta = nn.beta;
phi = 2 * pi / lambda * delta * z;
b = 2 * pi / lambda * beta * z;

fprintf('\nthickness z: %.1f mm',z*1000)
fprintf('\nenergy E = %.1f keV',energy/1000)
fprintf('\n%10s %10s','B','phase')
fprintf('\n%10f %10f',b,phi)

fprintf('\n')

  