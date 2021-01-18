function alpha = FineStructureConstant
%  The fine-structure constant, also known as Sommerfeld's constant,
%  commonly denoted by α (the Greek letter alpha). A fundamental physical
%  constant which quantifies the strength of the electromagnetic
%  interaction between elementary charged particles. It is a dimensionless
%  quantity related to the elementary charge e, which denotes the strength
%  of the coupling of an elementary charged particle with the
%  electromagnetic field, by the formula 4πε0ħcα = e2.

alpha = ElementaryCharge^2 / (4 * pi * VacuumPermittivity('F_m') * PlanckConstantReduced('Js') * SpeedOfLight );
% A^2 s^2 / ( F m^1 * J s * m s^-1 )
% A^2 s^2 / ( s^4 A^2 m^-2 kg^-1 m^1 * kg * m^2 s^2 s * m s^-1 )


