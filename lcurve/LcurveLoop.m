function  LcurveLoop(pts)
% L-curve method to determine regularization parameter for ESRF data.

if nargin<1, pts =3;end;

% Reco.m arguments.
% ESRF ID19 experimental setup parameters.
lambda = EnergyConverter(17.68);
distance = 1;
pixelsize = 0.69e-6;
padvalue = 0;
iterations = 0;
compute_correction = 0;

padding  = 0 ;  
alphamin = 2;
alphamax = 6;
dim1     = 2048;
dim2     = 2048;
ycrop    = 128;

% Folder to data.
ParentFolder = '/mnt/tomoraid3/tomo/ESRF_20100411_InhouseExperiment/';
ScanName     = {'CT_carbonMeshesPApet'; 'CT_graphiteMinePApet'; 'CT_meshesPApet'; 'CT_meshesPApet_B'};

dat = zeros(dim1,dim2);    

% Fourier space cooridnates for computation of the residual.
if padding >0,
    dimx        = padding*2^nextpow2(dim1);
    dimy        = padding*2^nextpow2(dim2);
    dat        = padarray(dat,[(dimx-dim1)/2,(dimy-dim2)/2],padvalue,'both');
else
    dimx    = dim1;
    dimy    = dim2;
end;
dimy = dimy - 2*ycrop;
[xi,eta]    = meshgrid(-1/2:1/(dimy):1/2-1/(dimy),-1/2:1/(dimx):1/2-1/(dimx));
xi          = fftshift(xi);
eta         = fftshift(eta);
lap         = xi.^2 + eta.^2;

% Loop over data.
for jj = 1:6, % Loop over distances within scan.
    for ii = 1:4, % Loop over samples defined in ScanName.
        % Get names of data files.
        tic;
        Prefix     = [char(ScanName(ii)) '_' num2str(jj) '_'];
        DataFolder = [ParentFolder char(ScanName(ii)) '/' Prefix '/'];
        RecoFolder = [DataFolder(1:end-1) 'PhaseRetrieval/'];
        files      = dir([DataFolder '/' char(ScanName(ii)) '_' num2str(jj) '_*.edf']);
        num_files  = numel(files);
        % Read flat and dark fields.
        dark       = edfread([DataFolder '/darkend0000.edf']);
        flat       = edfread([DataFolder '/refHST0000.edf']);
        %        flat(1,:,:) = edfread([DataFolder '/refHST0000.edf']);
        %        flat(2,:,:) = edfread([DataFolder '/refHST0500.edf']);
        %        flat(3,:,:) = edfread([DataFolder '/refHST1000.edf']);
        %        flat(4,:,:) = edfread([DataFolder '/refHST1500.edf']);        
        % Read data.
        dat   = edfread([DataFolder char(files(1).name)]);
        im    = (dat - dark)./(flat - dark);
        % Crop image.
        dat = dat(:,1+ycrop:end-ycrop);
        im  = im(:,1+ycrop:end-ycrop);
        tread = toc;
        fprintf(1,['Read first of ' num2str(num_files) ' in ' DataFolder '\n']);
        % Reconstruct phase (NLO) for different alpha.
        tic;
        for ii=1:pts;
            alpha(ii) = alphamin + (alphamax - alphamin)*(ii-1)/(pts-1);
            phi = Reco(im,alpha(ii),lambda,distance,pixelsize,padding+ ...
                       1,padvalue,iterations,compute_correction);
            phi = phi(:,:,1);
            phi_alpha(:,:,ii)  = phi;
            % Norm of the solution.
            sol_norm1(ii) = log10(norm(phi));
            % Residual norm of the solution.
            res_norm1(ii) = (norm(ifft2(lap.*fft2(phi))-dat));
        end;
        trec = toc;
        % Lcurve method.
        tic;
        
        % Compute curvature.
        sn1_d1 = diff(abs(sol_norm1),1);
        rn1_d1 = diff(abs(res_norm1),1);
        sn1_d2 = diff(abs(sol_norm1),2);
        rn1_d2 = diff(abs(res_norm1),2);
        sn1_d1 = interp1(1.5:pts,sn1_d1,2:pts-1);
        rn1_d1 = interp1(1.5:pts,rn1_d1,2:pts-1);
        k      = abs((rn1_d1.*sn1_d2 - rn1_d2.*sn1_d1)./(sn1_d1.^2 + rn1_d1.^2).^1.5);
        i_kmax = find(k==max(k));
        tcur   = toc;
        % Print parameters and save data.
        fprintf(1,['Reconstructed ' num2str(pts) ' phase maps for alpha in [' ...
                   num2str(alphamin) ',' num2str(alphamax) '] ']);
        fprintf(['(Elapsed time: %gs (reading), %gs (reconstruction). ' ...
                 '%gs (L-curve).\n'],tread,trec,tcur);
        fprintf(['At index %i where alpha = %g the curvature is maximum ' ...
         'with k = %g\n'],i_kmax,alpha(i_kmax),k(i_kmax));
        % Write binary stack of phi(x,y,alpha).
        unix(['mkdir ' RecoFolder ' >&/dev/null']);
        mexVolWrite([RecoFolder 'PhaseLO_vs_alpha_0000.vol'],phi_alpha,'float32')
        % Lcurve: Leading order.
        h = figure('visible','off');
        plot(res_norm1,sol_norm1,'.-');
        for kk=1:pts
            text(res_norm1(kk),sol_norm1(kk),['  (' num2str(kk) ', ' ...
                    num2str(alpha(kk)) ', ' num2str(10^-alpha(kk)) ')']);
        end
        saveas(h,[RecoFolder 'Lcurve_LO.jpg'],'jpg');
        % Curvature of L-curve (LO) vs alpha
        h = figure('visible','off');
        plot(2:pts-1,k,'.-');
        saveas(h,[RecoFolder 'LcurveCurvature_LO.jpg'],'jpg');
        % Phase at optimum.
        imwrite(normat(phi_alpha(:,:,i_kmax)),[RecoFolder ...
                            'PhaseAtOptimalAlpha.jpg'],'jpg','Quality',100);
        fid = fopen([RecoFolder 'AlphaOptimum.txt'],'w');
        fprintf(fid,num2str(alpha(i_kmax)));
        fclose(fid);
    end;
end;
    

