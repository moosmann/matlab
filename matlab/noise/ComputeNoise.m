format long
fprintf('\nPoisson statistics\n')
if 0
    % Read refs.
    ref0 = readstack(pwd,'ref*_0000','edf');
    % Define ROI.
    x = 500:2048-500;
    y = x;
    ref = ref0(x,y,:);
    % Normalize refs.
    for ii = 15:-1:1
        m = ref(:,:,ii);
        mm(ii) = mean(m(:));
        fprintf('%7omo3.6f\n',mm(ii));
    end
end
% Mean and var of refs.
refmean = mean(ref,3);
refvar = var(ref,1,3);
% Mean and var refs.
refnmean = mean(refn,3);
refnvar = var(refn,1,3);
% Means of means and vars.
meanrefmean = mean(refmean(:));
meanrefvar  = mean(refvar(:));
% Multiplication factor.
mp  = meanrefvar/meanrefmean;
% Print.
fprintf('2D mean of 1D mean along 3rd dim: %.5f\n',meanrefmean)
fprintf('2D mean of 1D var along 3rd dim:  %.5f\n',meanrefvar)
fprintf('Multiplication factor: %.5f\n',mp)
fprintf('Poisson mean of refs: %.5f\n',meanrefmean/mp)
