
% Phase retrieval.
[lo,nlo,phase,int,ldp,alphi]=Schwing(150,20,1);
%Good values for [sorder,smaxfac]=[150,20]. Then the differences between
%mass regularization and Schwinger regularization are neglectable.
% Phase maps.
if 0,
    ishow(nlo);
    ishow(alphi(:,:,2));end;
% Difference map.
ishow(abs(nlo-alphi(:,:,2)));
% Line plots.
if 1,
    LinePlots(phase,cat(3,lo,nlo));
    LinePlots(phase,alphi);end;
% Mean errors between retrieved and exact phase.
me(phase,lo);
me(phase,lo+nlo);
me(phase,alphi(:,:,1));
me(phase,alphi(:,:,1)+alphi(:,:,2));
% Mean error between Schwinger and mass regularization.
me(lo,alphi(:,:,1));
me(nlo,alphi(:,:,2));
me(nlo+lo,alphi(:,:,2)+alphi(:,:,1));
