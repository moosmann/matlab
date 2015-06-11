function eV = kJpermol_to_eV(kJ_per_mol)
% Convert kJ/mol to eV.
% 1 kJ·mol−1 is equal to 0.239 kcal·mol−1 or 1.04×10-2 eV per particle.
% Für ein einzelnes Elektron wird die Ionisierungsenergie in eV/Atom angegeben, 
% für 1 Mol Elektronen aber in kJ/mol. Der Umrechnungsfaktor ergibt sich aus 
% der Umrechnung zwischen eV und kJ sowie der Avogadro-Konstante N_A zu:
% 1 eV/particle = 96.485307 kJ/mol. 

eV = kJ_per_mol / 96.485307;
