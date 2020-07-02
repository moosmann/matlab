

% Flux measurements, commissioning, start: Wednesday, 2018-11-07
% PETRA III: 40 Bunch mode, 90 mA top-up

%Frontend Filter: 300 µm CVD + 50 µm Cu__  
% x-ray transmission at 30000 eV: 0.59__

%% Jupyter notebook issues
% - conversion factor flux(ss)/diode_avg differs from energy to energy
% - same values at E = 20 keV for 2x2 and 4x4

% FABIAN: ampfactor

%% 30 keV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter: 300 µm CVD + 50 µm Cu

nn = 1;
s(nn).energy = 30000;
ampfactor = 10000.000;
ss = 0;

ss = ss + 1;
s(nn).slit(ss) = 1*1;% mm^2:
s(nn).flux(ss) = 525542981736.40088;

ss = ss + 1;
s(nn).slit(ss) = 2 * 2; % mm^2
s(nn).flux(ss) = 1940990239781.5906;

ss = ss + 1;
s(nn).slit(ss) = 3 * 3; % mm^2
s(nn).flux(ss) = 4011265338705.3613;

ss = ss + 1;
s(nn).slit(ss) = 4 * 4; % mm^2
s(nn).flux(ss) = 6468891040609.7197;

ss = ss + 1;
s(nn).slit(ss) = 5 * 5; % mm^2
s(nn).flux(ss) = 7197470294599.7373;

% last: open
ss = ss + 1;
s(nn).slit(ss) = 5 * 5; % mm^2
s(nn).flux(ss) = 7196073558201.9785;

%first: open
ss = ss + 1;
s(nn).slit(ss) = 5 * 5; % mm^2
s(nn).flux(ss) =  7195894765156.3301;


%% E 40 keV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = nn + 1;
s(nn).energy = 40000;
ss = 0;

s(nn).time = 3;

ss = ss + 1; 
s(nn).slit(ss) = 7 * 7;% mm^2
s(nn).flux(ss) = 10488633252830.309;

ss = ss + 1;
s(nn).slit(ss) = 6 * 6; % mm^2
s(nn).flux(ss) = 10488633252830.309;

ss = ss + 1;
s(nn).slit(ss) = 5 * 5; % mm^2
s(nn).flux(ss) = 7599308484093.8027;

ss = ss + 1; 
s(nn).slit(ss) = 4 * 4; % mm^2
s(nn).flux(ss) = 5272609387773.2314;

ss = ss + 1; 
s(nn).slit(ss) = 3 * 3; % mm^2
s(nn).flux(ss) = 3145347014571.5469;

ss = ss + 1;
s(nn).slit(ss) = 2 * 2; % mm^2
s(nn).flux(ss) = 1469944935441.7734;

ss = ss + 1; 
s(nn).slit(ss) = 1 * 1; % mm^2
s(nn).flux(ss) = 374143014711.56433;

%% synchroload beamtime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = nn + 1;
s(nn).energy = 45000;
s(nn).flux = 1095734443447;
fd2 = 2.5655e+17 / 1000^3; %photons / s / m^2
s(nn).slit = 4.271000; % mm^2



%% Plot
figure( 'Name', 'Flux density vs slit' )
l = cell( [1, 3] );
for nn = 1:3
    flux_density = s(nn).flux  ./ (s(nn).slit * 1e-6);
    plot( s(nn).slit,  flux_density, '-x' )
    hold on
    l{nn} = sprintf( '%u keV', s(nn).energy / 1000 );
end
xlabel( 'slit area / mm^2' )
ylabel( 'flux density / ( photons / s / m^2' )
legend( l )


%% Plot
figure( 'Name', 'Flux vs slit' )
l = cell( [1, 3] );
for nn = 1:3
    flux = s(nn).flux  ;
    plot( s(nn).slit,  flux, '-x' )
    hold on
    l{nn} = sprintf( '%u keV', s(nn).energy / 1000 );
end
xlabel( 'slit area / mm^2' )
ylabel( 'flux  / ( photons / s )' )
legend( l )

%% E 50 keV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nn = nn + 1;
%energy(nn) = 50000; % eV
% Filter: 300 µm CVD + 50 µm Cu and 4mm GC
% x-ray transmission at 50000 eV: 0.73
