%Stadtwerke Tarifberechnungen


%%Tarif 1
t1 = 0.077;
g1 = 3.35;

%% Tarif 2
t2 = 0.054;
g2 = 13.50;

kwh1 = 4719;
kwh2 = 2886;

v = kwh2*t2;
g = 3*g2;
solln = v +g;
sollb = solln*1.19;
fprintf('\nTarif 2: %g + %g = %g (netto), %g (brutto)\n',v,g,solln,sollb)
v = kwh2*t1;
g = 3*g1;
solln = v +g;
sollb = solln*1.19;
fprintf('\nTarif 1: %g + %g = %g (netto), %g (brutto)\n',v,g,solln,sollb)


