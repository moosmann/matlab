% Script to simulate upcoming experiment at TopoTomo to get some
% preknowledge to optimize setup


PhotonCounts = 0;
compute_correction = 1;
blurring = 0;
alpha = 12;
EnergyDistancePixelsize = [25 0.01 0.36e-6];
compute_container = 1;

tic,
[phase,int,phiCTF,phiBRO,namestring]=Xenopus(PhotonCounts,compute_correction,blurring,alpha,EnergyDistancePixelsize);
   assignin('base',['xeno_pha'],phase);
   assignin('base',['xeno_int',namestring],int);
   assignin('base',['xeno_ctf',namestring],phiCTF);
   assignin('base',['xeno_bro',namestring],phiBRO);
t=toc;

if compute_container,
    compute_correction = 1;
    [phase,int,phiCTF,phiBRO]= Container(PhotonCounts,compute_correction, ...
                                         alpha,EnergyDistancePixelsize); 
   assignin('base',['cont_pha'],phase);
   assignin('base',['cont_int',namestring],int);
   phiCTF=phiCTF(:,:,1);assignin('base',['cont_ctf',namestring],phiCTF);
   phiBRO=phiBRO(:,:,1);assignin('base',['cont_bro',namestring],phiBRO);
end;

fprintf(1,'Elapsed time for computation: %gs (%s)\n',t,namestring);
clear phase int phiCTF phiBRO namestring PhotonCounts compute_correction ...
    compute_container blurring alpha EnergyDistancePixelsize t ans;  
 