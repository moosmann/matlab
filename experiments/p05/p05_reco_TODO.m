% TODO list for P05 reconstruction pipeline 'p05_reco'

% TODO: visual output for phase retrieval: phase map, fft(int), etc
% TODO: choose regulatrization parameter for phase retrieval via interactive mode 
% TODO: if visual output none, reduce correlation to chosen method
% TODO: physically consistent attenutation values of reconstructed slice 
% TODO: check raw roi, transpose, rot90 for KIT and EHD camera
% TODO: get rid of transpose when reading image files
% TODO: vertical stitching
% TODO: volume shape for excentric rot axis
% TODO: Interactive plot mode for wavelet-fft fing filter
% TODO: Improve FilterStripesCombinedWaveletFFT
% TODO: check proj-flat correlation measure for absolute values
% TODO: SSIM: include Gaussian blur filter, test
% TODO: large data set management: parloop, memory, etc
% TODO: stitching: optimize and refactor, memory efficency, interpolation method
% TODO: interactive loop over tomo slices for different phase retrieval parameter
% TODO: automatic determination of rot center
% TODO: output file format option: 8-bit, 16-bit for sino, phase, proj
% TODO: additional padding schemes for FBP filter
% TODO: read sinogram option
% TODO: set photometric tag for tif files w/o one, turn on respective warning
% TODO: GPU phase retrieval: parfor-loop requires memory managment
% TODO: median filter width of ring filter dependence on binning
% TODO: check offset: proj correlation for rotation axis determination
% TODO: check offset: flat/proj correlation 
% TODO: inverse Gaussian filter for phase retrieval, VZC theorem
% TODO: flat-flat correlation and averaging before proj-flat correlation
% TODO: proj-flat correlation: take negative logarithm or not?
% TODO: Check impact of symmetric FBP ifft option        
