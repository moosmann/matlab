function epsilon0 = VacuumPermittivity(units)
% The vacuum permittivity ε0, also called permittivity of free
% space or electric constant or the distributed capacitance of the vacuum,
% is an ideal, (baseline) physical constant, which is the value of the
% absolute dielectric permittivity of classical vacuum.
%
% ARGUMENTS
% units: string, 'F_m': Farady per metre or 'e2_eVm': electron rest mass
%   squared per electron Volt per metre or 'e2_GeVfm': electron rest mass
%   squared per Giga electron Volt per femtometre.

if nargin < 1
    units = 'F_m';
end

switch units
    case 'F_m'
        epsilon0 = 8.8541878128e-12; % F⋅m−1;
    case 'e2_eVm'
        epsilon0 = 55.26349406e-6; % e2⋅eV−1⋅m−1
    case 'e2_GeVfm'
        epsilon0 = 55.26349406; % e2⋅GeV−1⋅fm−1
end