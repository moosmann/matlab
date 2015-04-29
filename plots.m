function plots(stack,slice,padding,alpha,beta,fig);
%plot phase reconstruction

%compute Bronnikov and corrections
[phi0,phi1,phi2] = rec(stack(:,:,slice),padding,alpha,beta);

%figures
if fig==1
%figure,imshow(proj_bins,[],'InitialMagnification','fit'),colorbar;
figure('Name','Bronnkiov reconstruction: phi0'), ... 
    imshow(phi0,[],'InitialMagnification','fit'),colorbar;
end;

if fig==2
figure('Name','Correction: phi1'),  ...
    imshow(phi1+phi2,[],'InitialMagnification','fit'),colorbar;

end;

if fig==3
figure('Name','Intensity, Bronnikov, All corrections'), ...
  imshow([phi0,phi1+phi2],[],'InitialMagnification','fit');
end;


if fig==4
figure('Name','Intensity, Bronnikov, All corrections'), ...
  imshow([normat(dataslice),normat(phi0),normat(phi1+phi2)],[],'InitialMagnification','fit');
end;

