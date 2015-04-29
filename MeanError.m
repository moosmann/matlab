function [m0,m1] = MeanError(phase_exact,phase_rec_0,phase_rec_1,phase_rec_2)

m0 = mean(mean(abs(normat(phase_exact)-normat(phase_rec_0))));
m1 = mean(mean(abs(normat(phase_exact)-normat(phase_rec_0+phase_rec_1+phase_rec_2))));
