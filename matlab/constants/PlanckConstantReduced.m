function hbar = PlanckConstantReduced(units)
% Return the reduced Planck constant in electron volts and seconds or in
% Joule and seconds if 'units' is 'eVs' or 'Js', respectively.

if nargin < 1
    units = 'eVs';
end

hbar = PlanckConstant(units) / 2 / pi;

