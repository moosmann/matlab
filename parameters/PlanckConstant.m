function h = PlanckConstant(units)
% Return the Planck constant in electron volts and seconds or in Joule and
% seconds if 'units' is 'eVs' or 'Js', respectively.

if nargin < 1
    units = 'eVs';
end


switch units
    case 'eVs'
        h = 4.135667662e-15; % eV * s
    case 'Js'
        h = 6.626070040e-34;  % J * s
end

