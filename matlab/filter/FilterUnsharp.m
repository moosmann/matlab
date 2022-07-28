function [v, p] = FilterUnsharp(im_or_folder, radius, amount, threshold, outpath)
% Sharpen image using the unsharpen mask.
%
% im_or_folder: 2D image or folder name (string). If input is an image,
%   there will be an output. If input is a string, all tiff images will be
%   filter and saved in the folder [im_or_folder '_sharpenR*A*T*].
% radius: Standard deviation of Gaussian lowpass filter, 1 (default) |
%   positive number 
% amount: Strength of sharpening effect: 0.8 (default) | number
% threshold: Minimum contrast required for a pixel to be considered an edge
%   pixel 0 (default) | number in the range [0, 1]
%

if nargin < 2
    radius = 1;
end
if nargin < 3
    amount = 0.8;
end
if nargin < 4
    threshold = 0;
end
if nargin < 5
    outpath = '/beegfs/desy/user/moosmanj/filter_unsharp';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if image or folder sting is provided
if ischar(im_or_folder)
    if nargout > 0
        error( 'There will be output if a folder is given as first argument')
    end
    d = dir( [im_or_folder filesep '*tif']);
else
    im = im_or_folder;
    d = 1;
    nn = 1;
    num_par = numel(radius)*numel(amount)*numel(threshold);    
    v = zeros( [size(im,1), size(im,2), num_par], 'single');
    p(num_par).radius = [];
    p(num_par).amount = [];
    p(num_par).threshold = [];
end

% Loop over images
for s = 1:numel(d)
    if numel(d) > 1
        fn = sprintf( '%s/%s', d(s).folder, d(s).name);
        im = imread(fn);
    end
    % Loops over parameters: radiius, amount, threshold
    for r = 1:numel(radius)
        rr = radius(r);
        for a = 1:numel(amount)
            aa = amount(a);
            for t = 1:numel(threshold)
                tt = threshold(t);
                % Sharpen
                ims = imsharpen( im, 'Radius', rr, 'Amount', aa, 'Threshold', tt);
                % Save image first arguments is a folder
                if numel(d) > 1
                    % Output filename
                    fn2 = sprintf( '%s_sharpen/R%05.2f_A%04.2f_T%04.2f', d(s).folder, rr, aa, tt);
                    fn2 = regexprep(fn2,'\.','p');
                    CheckAndMakePath(fn2);
                    fn2 = sprintf( '%s/%s', fn2, d(s).name);
                    % Save image
                    write32bitTIF( ims, fn2);
                else
                    v(:,:,nn) = ims;
                    p(nn).radius = rr;
                    p(nn).amount = aa;
                    p(nn).threshold = tt;
                    fprintf( '\n%3u: %10f %10f %10f', nn, rr, aa, tt);
                    
                    fn2 = sprintf( '%s/im%04u_R%05.2f_A%04.2f_T%04.2f', outpath, nn, rr, aa, tt);
                    fn2 = regexprep(fn2,'\.','p');
                    fn2 = sprintf( '%s.tif', fn2);
                    % Save image
                    write32bitTIF(fn2, ims);
                    
                    nn = nn + 1;
                end
            end
        end
    end
end
