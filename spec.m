function spec(m)
    
    mf = fft(m,[],1);
    display(mf(1:10,256));
    lf = fft(m(:,256));
    display(lf(1:10));
    mf(1:10,256) == lf(1:10)
    