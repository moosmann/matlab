function lambda = ElectronComptonWavelength
% Compton wave length, lambda, of an electron:
%   lambda = h / m_e / c = 2.42631023867(73)×10^-12 m

lambda = PlanckConstant('Js') / ElectronRestMass('kg') / SpeedOfLight;

