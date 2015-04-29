function stack = PhaseSinoFilter(stack,PhaseMethod,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold)
% Combined sinogram filtering and phase retrieval
% Written by Julian Moosmann, last modified: 2013-09-18

% Fourier space filter for phase retrieval
phaseFilter = PhaseFilter(PhaseMethod,size(stack),EnergyDistancePixelsize,RegPar,BinaryFilterThreshold);
%% Sinogram filtering
stack = fft(fft(stack,[],2),[],3);
stack(:,:,1) = median(stack(:,:,[1:3 end-1:end]),3);
%% Phase retrieval
stack = real(ifft2(bsxfun(@times,fft(ifft(stack,[],3),[],1),phaseFilter)));

end
