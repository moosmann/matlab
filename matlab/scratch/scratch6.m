nn = 1;

            ind_range = (nn - 1 ) + (1:dpc_steps);
            ims = proj(:,:,ind_range);
            ft_ims = fft( ims, [], 3 );
            % Phase
            im = wrap( angle( ft_ims(:,:,2)) - angle( ft_flat_ims(:,:,2)) );
            dpc_phase(:,:,nn) = Binning( im, dpc_bin ) / dpc_bin^2 ;
            % Dark
            im = (abs( ft_ims(:,:,2) ) ./ abs( ft_ims(:,:,1) ) ) ./ (abs(ft_flat_ims(:,:,2))./abs(ft_flat_ims(:,:,1)));
            dpc_dark(:,:,nn) = Binning( im, dpc_bin ) / dpc_bin^2;
            % Attenuation
            im = abs( ft_ims(:,:,1)) ./ abs( ft_flat_ims(:,:,1) );
            dpc_att(:,:,nn) = Binning( im, dpc_bin ) / dpc_bin^2;