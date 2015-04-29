function stack = DataProc32ID_FilterFS(stack,StackString)
% Function to filter stack, called from main function DataProc32ID.

%% Start parallel Loop
parfor nn = 1:size(stack,1)% dim1 = size(stack,1);
    im = squeeze(stack(nn,:,:));
    im = fft2(im);
    im(2:end,1) = 0; im(1,2:end) = 0;
    stack(nn,:,:) = real(ifft2(im));
    %stack(nn,:,:) = real(ifft2(mask.*fft2(squeeze(stack(nn,:,:)))));
end
fprintf('FS Filtering: %s. ',StackString)
