function ComptonScattering()
% Compton Scattering analysis based on the Klein Nishina differential
% cross section.

fov = 0.003; %m
d = sqrt(2) * fov;
fprintf( '\nFOV diagonal : %f cm', d )


%energy = ([5 15 30 60 80 90 150 200] * 1e3)'; % eV
energy = ([ 60 80 90] * 1e3)'; % eV
%energy = ([2.75 60 511 1460 10000] * 1e3energy)'; % eV
theta = (0:1:100) / 100 * pi;
for n = numel( energy):-1:1
    l{n} = sprintf( '%g keV', energy(n) / 1000);
end

fprintf( '\n plot energies / keV: [' )
fprintf( '%u ', energy / 1000 )
fprintf( ']' )
fprintf( '\n plot angles: = [%u %u]', theta(1), theta(end) )

figure( 'Name', 'Ratio of photon energy after and before collision' )
t = theta * 180;
p = P( energy, theta );
plot( t / pi, p )
xlabel( 'angle / degrees' )
ylabel( 'energy / keV' )
legend( l )


figure( 'Name', 'Klein-Nishina differential cross section' )
t = theta * 180;
dcs = KN_diff_cs( energy, theta );
plot( t / pi, dcs )
xlabel( 'angle / degrees' )
ylabel( 'energy / keV' )
legend( l )

figure( 'Name', 'Klein-Nishina differential cross section: polarplot' )
t = [-theta(end:-1:1) theta];
dcs = KN_diff_cs( energy, t );
polarplot( t , dcs )
legend( l )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure( 'Name', 'Klein-Nishina forward scattering' )
z = [0.01 1.4]; %m
energy = [30:10:100] * 1e3; % eV
oa = atan( d ./ z );
fprintf( '\nz range: [%g %g]', z)
fprintf( '\nopening angle range: [%g %g]', oa)
a = oa(1):-range(oa)/100:oa(end);
zrange = d ./ tan( a);
r = zeros( [numel(energy), numel(a)]);
%fprintf( '\nenergy: %f keV', energy / 1000 )
for m = 1:numel(energy)
    em = energy( m );
    for n = 1:numel(a)
        an = a(n);
        [csp, cst] = KN_partial_cs( em, an, 10000 );
        r(m,n) = csp / cst;
    end
end
plot( zrange, r * 100 )
xlabel( 'z / m')
ylabel( 'partial / total cross-section in percent' )
title( 'forward scattering' )
for n = numel( energy):-1:1
    l{n} = sprintf( '%g keV', energy(n) / 1000);
end
legend( l )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry
fprintf( '\n' )
z = [0.01 1.4]; %m
fprintf( '\nz = [%f %f] cm', z*100)
oa = atan( d ./ z );
fprintf( '\nopening angle: [%f %f] degree', oa * 180 / pi )


fprintf( '\n' )
% Partial cross section
%[csp, cst] = KN_partial_cs( 30e3, 1, 1000 );

energy = 10e3;
opening_angle = 180;
[csp, cst] = KN_partial_cs( energy, opening_angle, 10000 );
fprintf( '\nopening angle: %f', opening_angle )
fprintf( '\nKN partial cross-section: %g', csp )
fprintf( '\nKN total cross-section: %g', cst )
fprintf( '\nKN ratio partial/total cross-section: %g', csp / cst )

cst2 = KN_total_cs( energy );
fprintf( '\nKN total cross-section: %g', cst2 )
fprintf( '\nKN ratio total/total cross-section: %g', cst2/cst )

fprintf( '\n' )



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ratio of photon energy before and after collision P = lambda / lambda'
function P = P(E_eV, theta)

P = 1 ./ ( 1 + ( E_eV / ElectronRestMass('eV_c2')) * ( 1 - cos(theta) ) );

% Klein-Nishina differntial cross section
function dcs = KN_diff_cs( E_eV, theta )

c = 1/2 * FineStructureConstant^2 * ElectronComptonWavelengthReduced^2;
p = P(E_eV, theta);
dcs = c * p.^2 .* ( p + p.^-1 - sin(theta).^2 );

% Klein-Nishina total cross section
function sigma = KN_total_cs( E_eV )
c = pi * FineStructureConstant^2 * ElectronComptonWavelengthReduced^2;
x = E_eV / ElectronRestMass('eV_c2');
sigma = c * x^-3 * ( 2*x*(2 + x*(1+x)*(8+x))/(1+2*x)^2 + ((x-2)*x-2)*log(1+2*x));

% Partial corss section of the forward propagating cone
function [cs_partial, cs_total] = KN_partial_cs( E_eV, opening_angle_degree, pts )
if nargin < 3
    pts = 10000;
end
oa = opening_angle_degree;
theta = (0:pts) / pts * pi;
dtheta = theta(2) - theta(1);
dcs = KN_diff_cs( E_eV, theta );


cs_total = 2 * pi * sum( dtheta * abs(cos(theta)) .* dcs, 2 );

m = theta <= opening_angle_degree * pi / 180;
cs_partial = 2 * pi * sum( dtheta * abs(cos(theta)) .* dcs .* m, 2 );



