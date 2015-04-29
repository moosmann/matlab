tieLO1 = RecoSchwingGauss(int,[16 sqrt(1024)*1024],EnergyConverter(30),0.5,0.6e-6,1,0,0);
tieLO2 = RecoSchwingGauss(int,[150 1024],EnergyConverter(30),0.5,0.6e-6,1,0,0);

LinePlots(phase,cat(3,tieLO1(:,:,1),tieLO2(:,:,1)),532,1:1000)