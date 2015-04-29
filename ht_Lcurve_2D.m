## Copyright (C) 2008 P. Cloetens
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

## ht_Lcurve_2D
## function [delta_beta_opt, lim1_opt] = ht_Lcurve_2D(me, re, delta_beta, lim1, resultdir)
##      L-curve method for several delta_beta values
##      optimum is choosen to correspond to the one yielding largest lim1_maxcurv value
##      To be tested

## Author: P. Cloetens <cloetens@esrf.fr>
##
## 2008-11-26 P. Cloetens <cloetens@esrf.fr>
## * Initial revision

function [delta_beta_opt, lim1_opt] = ht_Lcurve_2D(me, re, delta_beta, lim1, resultdir)
    lim1_maxcurv = lim1_minME = lim1_minspeed = zeros(size(delta_beta));
    for q = 1:length(delta_beta)
        [lim1_maxcurv(q), lim1_minME(q), lim1_minspeed(q), K(q,:), V(q,:)] = ht_Lcurve(me(q,:), re(q,:), lim1, \ 
                resultdir);
    endfor
    imagej(K); sleep(0.2); imagej_rename('K')
    imagej(V); sleep(0.2); imagej_rename('V')
    [lim1_opt, mapos] = max(lim1_maxcurv);
    delta_beta_opt = delta_beta(mapos);
    printf("Maximum value of lim1 reached for delta_beta value: %g\n", delta_beta_opt)
    printf("Corresponding lim1 value: %g\n", lim1_opt)
    
    figure;
    plot(delta_beta, lim1_maxcurv, 'rx-', delta_beta_opt, lim1_opt, 'gx')
    xlabel('delta_beta')
    ylabel('lim1_maxcurv')
    
    # Save results
    name = sprintf("%s/Lcurve_2D.mat", resultdir);
    save(name, 'delta_beta_opt', 'lim1_opt','K', 'V', 'lim1_maxcurv', 'lim1_minME', \
            'lim1_minspeed', 'me', 're', 'delta_beta', 'lim1')
endfunction
