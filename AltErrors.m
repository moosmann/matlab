npe = normat(pe);
npb = normat(p0);
npbc = normat(p0+p1+p2);
npe_mean = mean(npe(:));
npb_mean = mean(npb(:));
npbc_mean = mean(npbc(:));
npe = npe - npe_mean;
npb = npb - npb_mean;
npbc = npbc - npbc_mean;
errorpb = abs(npe-npb);
errorpbc = abs(npe-npbc);
mean_errorpb  = mean(errorpb(:));
mean_errorpbc = mean(errorpbc(:));
fprintf(1,['Mean real space error: BRO=%.8g, BROCOR=%g (%2.3f%%)\n'], ...
        mean_errorpb, mean_errorpbc,100*mean_errorpbc/mean_errorpb);