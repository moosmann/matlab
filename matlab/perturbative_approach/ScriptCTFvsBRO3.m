%Situation: Phasenhub ist nicht klein. Waehle phi_max=1.

%Gehe zu grossem lambda*z/pixelsize fuer besseren Intensitaetskontrast
%bzw. besseres Signalrauschverhaeltnis.
ca
%E=10,z=50: zero-crossing~12
if 1,
if 0,ishow(int_1p00_E10_z0p50_c00000),
ishow(int_1p00_E10_z0p50_c20000),
ishow(int_1p00_E10_z0p50_c01000),
ishow(ctf_1p00_E10_z0p50_c00000,1),
end;
ishow(ctf_1p00_E30_z0p10_c00000,1)
ishow(ctf_1p00_E30_z0p10_c08000,1)

ishow(ctf_1p00_E10_z0p50_c00000,1)
ishow(ctf_1p00_E10_z0p50_c08000,1)

ishow(bro_1p00_E30_z0p10_c00000,1)
ishow(bro_1p00_E30_z0p10_c08000,1)

ishow(bro_1p00_E10_z0p50_c00000,1),
ishow(bro_1p00_E10_z0p50_c08000,1)

end;
