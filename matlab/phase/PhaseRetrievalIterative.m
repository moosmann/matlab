%flatcor_path = '/asap3/petra3/gpfs/p07/2019/data/11006991/processed/syn008_Ti_12w_47R_z0030_t135ms/flat_corrected/rawBin4';

flatcor_path = '/asap3/petra3/gpfs/p05/2020/data/11009667/processed/embl_056_200924_dist_1400_zshift_0p2_010/flat_corrected/rawBin2/';
phase_map_path = '/asap3/petra3/gpfs/p05/2020/data/11009667/processed/embl_056_200924_dist_1400_zshift_0p2_010/phase_map/rawBin2/tie_regPar1p10';
name = 'embl_056_200924_dist_1400_zshift_0p2_010';
num_iter = 50;
constraint = 0;
dtype = 'double';
show_fig = 0;


ds_att = dir( [flatcor_path filesep '*.tif'] );
ds_pha = dir( [phase_map_path filesep '*.tif'] );
num_im = numel( ds_att );
num_im_half = round(num_im/2);

outpath = '/asap3/petra3/gpfs/p05/2020/data/11009667/processed/embl_056_200924_dist_1400_zshift_0p2_010/phase_map_iterative';
CheckAndMakePath( outpath )

for p = 1:num_im
    fprintf( ' %u', p )
    %for p = num_im_half
    
    % Intensity / attenuation
    filename = sprintf( '%s/%s', flatcor_path,  ds_att(p).name );
    int = imread( filename );
    
    % phase map
    if numel( ds_pha ) ~= num_im
        error('phase maps and flat-corrected projecitons do not match');
    end
    filename = sprintf( '%s/%s', phase_map_path,  ds_pha(p).name );
    pha0 = imread( filename );
    pha0 = -pha0;
    
    % Type conversion
    switch dtype
        case 'single'
            pha0 = single( pha0 );
            int = single( int );
        case 'double'
            pha0 = double( pha0 );
            int = double( int );
    end
    if isempty( dtype )
        dtype = class( int );
    end
    
    % Display flatcor and tie phase map
    if show_fig
        im_size = size( int );
        if exist( 'h1' , 'var' ) && isvalid( h1 )
            figure(h1)
        else
            h1 = figure( 'Name', 'flatcor and phase map' );
        end
        subplot(2,1,1)
        imsc( int )
        title(sprintf('int #%u', p))
        colorbar
        axis equal tight
        
        subplot(2,1,2)
        imsc( pha0 )
        title(sprintf('pha non iterative #%u', p))
        colorbar
        axis equal tight
        
        drawnow
    end
    
    % Parameters.
    energy      = 30000; % keV
    lambda = E_to_lambda(energy); % m
    distance    = 1.399999500000000; % m
    pixelsize   =  2.1400001e-06; % m
    fresnel_number = pixelsize^2 / lambda / distance;
    if show_fig
        fprintf( 'Parameter')
        fprintf( '\n energy : %f keV', energy / 1000 )
        fprintf( '\n z : %f m', distance )
        fprintf( '\n pixe size : %f micron', pixelsize * 1e6 )
        fprintf( '\n Fresnel number : lambda * z / dx^2 =  %f', fresnel_number)
    end
    
    % Fourier cooridnates 1D
    normalize = 1;
    xi1 = FrequencyVector(2*im_size(2), dtype, normalize);
    eta1 = FrequencyVector(2*im_size(1), dtype, normalize);
    [xi,eta] = meshgrid(xi1, eta1);
    if show_fig
        fprintf( '\n Fourier cooridantes :[%f, %f, ..., %f, %f; %f, %f, ... , %f, %f] ', xi1(1:2), xi1(im_size(2) + (-1:2)), xi1(end-1:end) )
        fprintf( '\n Fourier cooridantes :[%f, %f, ..., %f, %f; %f, %f, ... , %f, %f] ', eta1(1:2), eta1(im_size(1) + (-1:2)), eta1(end-1:end) )
    end
    
    % Propagator
    fp_cpu = exp( -1i * pi / 2 / fresnel_number * (xi.^2 + eta.^2) );
    fp = gpuArray( fp_cpu );
    
    % Figure and plot range
    x = 1:im_size(1);
    y = 1:im_size(2);
    if show_fig
        if exist( 'h2' , 'var' ) && isvalid( h2 )
            figure(h2)
        else
            h2 = figure( 'Name', 'iterative phase retrieval' );
        end
    end
    
    %% phase only
    if 1
        
        % Init
        int_sqrt = gpuArray( padarray( sqrt( int ), im_size, 'symmetric', 'post') );
        int_sqrt_f = FilterOutlier( int_sqrt, [0.2 0.2], 'none' );
        int_sqrt_f_mean = mean2( int_sqrt_f );
        int_sqrt = int_sqrt / int_sqrt_f_mean;
        
        psic = padarray( exp( 1i * 0.8 * pha0 ), im_size, 'symmetric' , 'post');
        psi = gpuArray( psic );
        
        psi = gpuArray( complex( zeros( size( int_sqrt ), dtype ) ) );
        
        m = gpuArray( zeros( size( int_sqrt ), 'logical' ) );
        
        if show_fig
        fprintf( '\n\nStart iteration' )
        fprintf( '\n' )
        end
        tic
        for num = 1:num_iter
            
            % Fresnel projector
            % forward propagation
            psi = int_sqrt .* exp( 1i * angle( ifft2( fp .* fft2( psi ) ) ) );
            % backward propagation
            psi = ifft2( 1 ./ fp .* fft2( psi ) );
            
            % phase and amplitude
            pha = angle( psi );
            amp = abs( psi );
            
            if show_fig && ( mod(num, round(num_iter / 10) ) == 0 )
                fprintf( '\n' )
                domain( pha, 1, 'phase before' )
                domain( pha, 1, 'amp before' )
            end
            
            % Constraint projector
            switch constraint
                case +1
                    m = pha < 0;
                    pha(m) = 0;
                    %psi = amp .* exp( 1i * pha );
                    psi = exp( 1i * pha );
                case -1
                    m = pha > 0;
                    pha(m) = 0;
                    %psi = amp * exp( 1i * pha );
                    psi = exp( 1i * pha );
                case 0
                    psi = exp( 1i * pha );
            end
            
            if show_fig && (mod(num, round(num_iter / 10) ) == 0)
                domain( pha, 1, 'phase after' )
                domain( pha, 1, 'amp after' )
            end
            
            
            % Plot
            if show_fig && ( mod(num, round(num_iter / 10) ) == 1 || num == num_iter )
                fprintf( '\n' )
                fprintf( ' %u', num )
                
                m = 2;
                n = 1;
                subplot(m,n,1)
                imsc( gather( pha(x,y) ) )
                title( sprintf( 'phase 0, iteration %u', num ) )
                colorbar
                axis equal tight
                
                subplot(m,n,2)
                imsc( gather( amp(x,y) ) )
                title( sprintf( 'amp 0, iteration %u', num ) )
                colorbar
                axis equal tight
                
                %             subplot(m,n,3)
                %             imsc( pha0 )
                %             title( sprintf( 'pha0' )  )
                %             colorbar
                %             axis equal tight
                
                drawnow
                pause( 0.5 )
            end
        end
    end
    
    %% Attenuation and phase retrieval
    if 0
        % Propagator
        fp0 = exp( -1i * pi * fresnel_number / 2 * (xi.^2 + eta.^2) );
        fp = gpuArray( fp0 );
        bp = gpuArray( 1 ./ fp0 );
        
        % Init
        attz_meas = gpuArray( padarray(int, im_size, 'symmetric', 'post') );
        phaz = gpuArray( zeros( size( attz_meas ), dtype ) );
        
        for num = 1:num_iter
            fprintf( ' %u', num )
            
            % Backpropagation
            phi = ifft2( bp .* fft2( exp( - attz_meas + 1i * phaz) ) ) ;
            att0 = -log( abs( phi ) );
            pha0 = angle( phi );
            
            
            % Constraints
            
            
            % Forward propagation
            phi = ifft2( fp .* fft2( exp( - att0 + 1i * pha0) ) ) ;
            attz = -log( abs( phi ) );
            phaz = angle( phi );
            
            % Constraints
            
            %     domain( att0 )
            %     domain( pha0 )
            %     domain( attz )
            %     domain( phaz )
            
            % Plot
            if mod(num, round(num_iter / 10) ) == 0
                fprintf( '\n' )
                
                m = 3; n = 2;
                subplot(m,n,1)
                imsc( gather( att0(x,y) ) )
                title( 'att0' )
                colorbar
                axis equal tight
                
                subplot(m,n,2)
                imsc( gather( pha0(x,y) ) )
                title( 'pha0' )
                colorbar
                axis equal tight
                
                subplot(m,n,3)
                imsc( gather( attz(x,y) ) )
                title( 'attz' )
                colorbar
                axis equal tight
                
                subplot(m,n,4)
                imsc( gather( phaz(x,y) ) )
                title( 'phaz' )
                colorbar
                axis equal tight
                
                subplot(m,n,5)
                imsc( gather( attz_meas(x,y) - attz(x,y)) )
                title( 'attz meas - attz' )
                colorbar
                axis equal tight
                
                subplot(m,n,6)
                imsc( gather( phaz(x,y) - pha0(x,y) ) )
                title( 'phaz - pha0' )
                colorbar
                axis equal tight
                
                drawnow
                %         pause(0.2)
            end
        end
    end
    phaz = gather( angle( psi(x,y) ) );
    
    % Save
    filename = sprintf( '%s/%s_%06u.tif', outpath, name, p );
    write32bitTIF( filename, phaz )
    
    if show_fig
        if exist( 'h3' , 'var' ) && isvalid( h3 )
            figure(h3)
        else
            h3 = figure( 'Name', 'phase maps' );
        end
        m = 2; n = 1;
        subplot(m,n,1)
        imsc( phaz )
        title( 'pha iterative' )
        colorbar
        axis equal tight
        
        subplot(m,n,2)
        imsc( pha0 )
        title( 'pha non iterative' )
        colorbar
        axis equal tight
        
        drawnow
    end
end
fprintf('\n')
toc
