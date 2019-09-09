function m = electron_rest_mass(units)
% Return the electron rest mass in kilogram ('kg'), in Joule over c^2
% ('J_c2'), in electron volts over c^2 ('eV_c2'), or in mega electron volts
% over c^2 ('MeV_c2').

if nargin < 1
    units = 'kg';
end

switch units
    case 'kg'
        m = 9.10938356e-31; % kg
    case 'J_c2'
        m = 8.18710565e-14; % J / c^2
    case 'eV_c2'
        m = 0.5109989461e6; % eV / c^2        
    case 'MeV_c2'
        m = 0.5109989461; % MeV / c^2        
end
