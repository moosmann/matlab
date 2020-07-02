function out = dose_in_vivo(fluxDensity_ph_per_s_mm2,energy_keV,absorptionLength_mm,scanTime_s,objectLength_mm,objectDensity_g_per_ccm,containerLength_mm)

if nargin < 1
    % 2-BM-B flux @ 30keV DMM: Henke: 8.0e13 ph/s/mrad^2/0.1%BW 
    % flux density = 8e13*10 ph/s/1%BW / (35m)^2 /mrad^2
    % = 6.54e11 ph/s/1%BW/mm^2 ~ 10^12 ph/s/mm^2
    fluxDensity_ph_per_s_mm2 = 10^12; % photons / s / mm ^2
end
if nargin < 2
    energy_keV = 30 ; % keV
end
if nargin < 3
    % for water: l_att = 1 / (4*pi*beta/lambda) = 3.10 cm at 30 keV
    absorptionLength_mm = 30.9;
end
if nargin < 4
    scanTime_s = 20 ; % s
end
if nargin < 5
    objectLength_mm = 1.0 ; % mm
end
if nargin < 6 
    objectDensity_g_per_ccm  = 0.997 ; % g / cm^3 = 1000 kg / m^3
end
if nargin < 7
    containerLength_mm = 12 ; % mm
end
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Attenuation
out.containerTransmission = exp( - (containerLength_mm -objectLength_mm)/2 / absorptionLength_mm) ;
out.objectAbsorption = 1-exp( - objectLength_mm/ absorptionLength_mm) ;

% Number of photons
out.absorbedPhotons = out.containerTransmission * out.objectAbsorption * fluxDensity_ph_per_s_mm2 * scanTime_s  * (objectLength_mm)^2;

% Cubic mass of object
out.cubeMass_g = objectLength_mm^3 * objectDensity_g_per_ccm / 1000 ;

% Dose
eV_to_J = 1.602176565e-19 ;
out.dose_Gy = out.absorbedPhotons * energy_keV*10^3 * eV_to_J ...
    / ( out.cubeMass_g / 1000 );

