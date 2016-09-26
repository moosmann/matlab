function [m0,m1,m0r,m1r] = SpectralAnlaysis(phase_exact, phase_rec_0, ...
                                            phase_rec_1,phase_rec_2)
    % Compute mean error of the reconstructed phase maps (Bronnikov,
    % Bronnikov plus Corrections) in 2D Fourier space.
    
    
    phe_f = fft2(normat(phase_exact));
    ph0_f = fft2(normat(phase_rec_0));
    ph1_f = fft2(normat(phase_rec_0 + phase_rec_1 + phase_rec_2));
    domain(phe_f);
    domain(ph0_f);
    domain(ph1_f);
       
    m0 = mean(mean(abs(phe_f - ph0_f)));
    m1 = mean(mean(abs(phe_f - ph1_f)));
    
    m0r = mean(mean(abs((real(phe_f)) - (real(ph0_f)))));
    m1r = mean(mean(abs((real(phe_f)) - (real(ph1_f)))));

    m0i = mean(mean(abs((imag(phe_f)) - (imag(ph0_f)))));
    m1i = mean(mean(abs((imag(phe_f)) - (imag(ph1_f)))));

    fprintf(1,['mean spectral error (abs): %1.2g (Bro), %1.2g (BroCor) ' ...
               '(%2.3f%%)\nmean spectral error (real): %1.2g (Bro), %1.2g ' ...
               '(BroCor)(%2.3f%%)\nmean spectral error (imag): %1.2g (Bro), %1.2g ' ...
               '(BroCor)(%2.3f%%)\n'],m0,m1,100*m1/m0,m0r,m1r,100*m1r/ ...
            m0r,m0i,m1i,100*m1i/m0i);
    