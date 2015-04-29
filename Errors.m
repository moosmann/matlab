function [mean_error_bro,mean_error_brocor,stadev_error_bro,stadev_error_brocor,error_bro,error_brocor] ...
     = Errors(exact_phase,reco_phase,pow)
    
if (nargin<1),exact_phase=0;end;
if (nargin<3),pow=1;end;
    
% 2D error maps in real space.
error_bro     = exact_phase - reco_phase(:,:,1);
error_brocor  = error_bro   - reco_phase(:,:,2);

% Scalar mean absolute error of the error maps.
mean_error_bro       = mean(abs(error_bro   (:)).^pow).^(1/pow);
mean_error_brocor    = mean(abs(error_brocor(:)).^pow).^(1/pow);
% Scalar standard deviation of the mean absolute error.
stadev_error_bro     = sqrt(var(abs(error_bro   (:))));
stadev_error_brocor  = sqrt(var(abs(error_brocor(:))));
% Print mean error and standard deviation.
fprintf(1,['Real space errors: MeanBro=%g, MeanBroCor=%g (%2.3f%%), StaBro=%g, StaBroCor=%g\n'], ...
 mean_error_bro,mean_error_brocor,100*mean_error_brocor/mean_error_bro,stadev_error_bro,stadev_error_brocor);
