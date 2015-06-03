clear all

roi = [0 1];
%roi = [0.5 0.85];

%% Parameters
% Acceleration voltage = 200 keV
energy_keV = 200;
lambda_m = EnergyConverter(energy_keV);
%Defocus: 3.0 micrometer under focus 
defocus_m = 3e-6; 
%Focal length: 2.7 mm
focalLength_m = 2.7e-3;
% Magnification: 25000
magnification = 2500;
%Detector pixel size: 16 micrometer
detectPixel_m = 16e-6;
effectPixel_m = detectPixel_m / magnification;
% Spherical aberration: 2.1 mm
cs_m = 2.1e-3;
% Chromatic aberration: 2.1 mm
cc_m = 2.2e-3; 
imSize = [210 250];
% Distance object lens
q_m = focalLength_m * ( magnification + 1 );

a = 41;

%lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
pfDefocus = pi * lambda_m * (a*defocus_m ) / effectPixel_m^2;
pfSpher = 0* pi * cs_m * lambda_m^3 / 2 / effectPixel_m^4;


%% Fourier coordinates
% 1D
outputPrecision = 'double';
xi  = FrequencyVector(imSize(2),outputPrecision,1);
eta = FrequencyVector(imSize(1),outputPrecision,1);
% 2D
[xi, eta]   = meshgrid(xi,eta);
% Function on 2D
xi = (xi.^2 + eta.^2);
xi = fftshift( xi );
arg = pfDefocus * xi - pfSpher * xi.^2;
%figure(1)
%imsc( 1./( 1 + abs( sin( arg ))))
%figure(2)

fprintf('\nDefocus prefactor: %g',pfDefocus)
fprintf('\nSpherical abberation prefactor: %g\n',pfSpher)
fprintf('Max argument: %g\n', max( arg(:) ) )

imagesc( sin(arg) )
colormap(gray)

%% Phase retrieval filter
regPar = 1.5;
pfctf = PhaseFilter('ctf', imSize, [energy_keV, a*defocus_m, effectPixel_m], regPar);
pftie = PhaseFilter('tie', imSize, [energy_keV, a*defocus_m, effectPixel_m], regPar);
pfqp = PhaseFilter('qp', imSize, [energy_keV, a*defocus_m, effectPixel_m], regPar, 0.15);
pfqp2 = PhaseFilter('qp', imSize, [energy_keV, a*defocus_m, effectPixel_m], regPar, 0.05);

%% Simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parentPath = '~/data/tem/JulianET/Phase retrieval';
pha = ReadMRC(sprintf('%s/balls_phantom.mrc', parentPath));
folderDose = 'Data set 1. Dose series';
folderDefocus = 'Data set 2. Defocus series';
dataFolder = folderDose;

% noise-free intensity
filename = sprintf('%s/%s/TEMimage_nonoise.mrc', parentPath, folderDose);
[int0, s, mi, ma, av] = ReadMRC(filename);
% function for plotting of Fourier space modulation
func = @(m) fftshift( log( 1+abs( fft2( SubtractMean(m) ) ) ) );

% loop over images
dataDir = dir( sprintf( '%s/%s/*.mrc', parentPath, dataFolder) );
for nn = numel(dataDir):-1:1
    str = dataDir(nn).name;
    digIndex = regexp(str,'\d');
    numList(nn) = str2double(str(digIndex));
end
% Remove NaNs
numList(isnan(numList)) = [];
numList = sort(numList);

figure(1)
%set(gca, 'XTickLabel', [],'XTick',[])
%tightfig;
x = round(1 + roi(1)*imSize(1)):round(roi(2)*imSize(1));
y = round(1 + roi(1)*imSize(2)):round(roi(2)*imSize(2));
imsc = @(map) imsc(map(x,y));

outPath = sprintf('%s/images/', parentPath);
%% Loop through all images sorted by digits, but exlcudes nonoise image
if 1
    for nn = 1:numel(numList)+1
        
        if nn < numel(numList)
            dataname = numList(nn);
            nnStr = sprintf('%06u', numList(nn));
            filename = sprintf( '%s/%s/TEMimage_%u.mrc', parentPath, dataFolder, dataname);
        else
            nnStr = 'nonoise';
            filename = sprintf( '%s/%s/TEMimage_%s.mrc', parentPath, dataFolder, nnStr);
        end
        
        
        [map, s] = ReadMRC(filename );
        
        subplot(3,2,1)
        imsc(map)
        title(dataname,'Interpreter', 'none')
        axis off tight image 
        imwrite(normat(map), sprintf('%s/int/int_%s.png', outPath, nnStr), 'png')
        
        subplot(3,2,2)
        imsc(func(map))
        title('log(1+abs(fft2(int)))')        
        axis off tight image        
        
        subplot(3,2,3)
        pha = ifft2( pftie.* fft2( map ) ); 
        imsc( pha )
        title('tie')
        axis off image
        imwrite(normat(pha), sprintf('%s/tie/tie_%s.png', outPath, nnStr), 'png')
        
        subplot(3,2,4)
        pha = ifft2( pfctf.* fft2( map ) );
        imsc( pha )        
        title('ctf')
        axis off image
        imwrite(normat(pha), sprintf('%s/ctf/ctf_%s.png', outPath, nnStr), 'png')
        
        subplot(3,2,5)
        pha = ifft2( pfqp.* fft2( map ) );
        imsc( pha )        
        title('qp')
        axis off image
        imwrite(normat(pha), sprintf('%s/qp/qp_%s.png', outPath, nnStr), 'png')
                
        subplot(3,2,6)
        pha = ifft2( pfqp2.* fft2( map ) );
        imsc( pha )        
        title('qp')
        axis off image
        imwrite(normat(pha), sprintf('%s/qp_2/qp_%s.png', outPath, nnStr), 'png')
        
        %pause(1)
        %tightfig;
        input('Hit Enter to continue')
    end    
end


%% Loop through all images in unsorted manner
if 0
    for nn=1:numel(dataDir)
        dataname = dataDir(nn).name;
        filename = sprintf( '%s/%s/%s', parentPath, dataFolder, dataname);
        [map, s] = ReadMRC(filename );                    
        
        subplot(3,2,1)
        imsc(map)
        title(dataname,'Interpreter', 'none')
        axis off tight image 
        
        subplot(3,2,2)
        imsc(func(map))
        title('log(1+abs(fft2(int)))')
        axis off tight image
        
        subplot(3,2,3)
        pha = ifft2( pftie.* fft2( map ) );
        imsc( pha )        
        title('tie')
        axis off image
        
        subplot(3,2,4)
        pha = ifft2( pfctf.* fft2( map ) );
        imsc( pha )        
        title('ctf')
        axis off image
        
        subplot(3,2,5)
        pha = ifft2( pfqp.* fft2( map ) );
        imsc( pha )        
        title('qp')
        axis off image
        
        subplot(3,2,6)
        pha = ifft2( pfqp2.* fft2( map ) );
        imsc( pha )
        title('qp')
        axis off image
        
        %pause(1)
        %tightfig;
        input('Hit Enter to continue')
    end
end
