%% Copyright (C) 2007 M. Langer
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

%% ht_Lcurve
%% function [lim1_maxcurv, lim1_minME, lim1_minspeed, K, V] = ht_Lcurve(me, re, lim1, resultdir, Lcurve_method, end_cond)
%%      L-curve method
%%      Find point of maximum curvature
%%      Find point of minimum speed
%%      Display log10 of regularization error versus model error and spline fit
%%      Display cuvature and speed as function of parameter lim1
%%      Create Lcurve.mat file with relevant output
%%
%%      outputs
%%      output 1: parameter value corresponding to maximum curvature
%%      output 2: parameter value corresponding to minimum model error
%%      output 3: parameter value corresponding to minimum speed
%%      output 4: oversampled curvature
%%      output 5: oversampled speed
%%      Relevant output is tored in Lcurve.mat
%%
%%      arguments
%%      argument 1: model error ( default: read from Lcurve.mat )
%%      argument 2: reglarization error ( default: read from Lcurve.mat )
%%      argument 3: parameter values ( default: read from Lcurve.mat )
%%      argument 4: resultdir ( default: current directory )
%%      argument 5: Lcurve_method: spline | polynomial ( default: 'spline' )
%%                  polynomial is not a good idea !!!
%%      argument 6: end condition for cubic spline, see function csape
%%                  ( default: 'not-a-knot' )
%%
%%      examples:
%%          ht_Lcurve; -> repeat L-curve method on existing results stored in Lcurve.mat
%%          ht_Lcurve(me, re, lim1); -> apply L-curve method on given arguments
%%          ht_Lcurve([], [], [], [], 'polynomial'); -> use polynomial fit instead of spline
%%                                                      not recommended
%%          ht_Lcurve([], [], [], [], [], 'variational'); -> use natural cubic spline
%%                                                           not recommended

%% Author: M. Langer <mlanger@esrf.fr>
%%
%% 2007-12-06 M. Langer <mlanger@esrf.fr> , P. Cloetens <cloetens@esrf.fr>
%% * Initial revision
%% 2008-04-03 PC
%% * clean up
%% * look for maximum curvature between values of lim1 that give min ME and min RE
%% * save results in Lcurve.mat
%% * plot curvature in separate window
%% * make function callable outside holotomo_slave
%% 2008-04-09 PC
%% * extended help
%% * calculate speed along curve and find point of minimum speed
%% * foresee polynomial fit next to spline
%% * foresee different end conditions for spline fit
%% * plot normalized curvature and speed in single figure
%% * change argument of fit and its subsampled version

function [lim1_maxcurv, lim1_minME, lim1_minspeed, K, V] = ht_Lcurve(me, re, lim1, resultdir, Lcurve_method, end_cond)

    if ~exist('me','var')
        me = [];
    end
    if ~exist('resultdir','var')
        resultdir = '';
    end
    if ~exist('Lcurve_method','var')
        Lcurve_method = '';
    end
    if ~exist('end_cond','var')
        end_cond = '';
    end
    
    if isempty(resultdir)
        resultdir = pwd;
    end
    if isempty(me)
        name = sprintf('%s/Lcurve.mat', resultdir);
        load(name);
    end
    if isempty(Lcurve_method)
        Lcurve_method = 'spline';
    end
    if isempty(end_cond)
        end_cond = 'not-a-knot';
    end
    
    disp('L-curve')

    mehat = log10(me);
    rehat = log10(re);
    N = 10;
    %t = 1:length(mehat);
    %ts = 1:1/N:length(mehat);
    t = linspace(0,1,length(mehat));
    ts = linspace(0,1,(length(mehat)-1)*N+1);
    limr = min(log10(lim1)):(max(log10(lim1))-min(log10(lim1)))/(length(lim1)-1)/N:max(log10(lim1));
    
    if strcmp(Lcurve_method,'spline')
        %corner by spline curvature
        D = [0 0 0 0; 3 0 0 0; 0 2 0 0; 0 0 1 0];
        
        px = csape(t,mehat,end_cond);
        py = csape(t,rehat,end_cond);
        
        [yy, P] = unmkpp(py);
        Y = ppval(py,ts);
        DP = P*D'; Dpy = mkpp(yy, DP); DY = ppval(Dpy,ts);
        D2P = DP*D'; D2py = mkpp(yy, D2P); D2Y = ppval(D2py,ts);
        
        [xx, P] = unmkpp(px);
        X = ppval(px,ts);
        DP = P*D'; Dpx = mkpp(xx, DP); DX = ppval(Dpx,ts);
        D2P = DP*D'; D2px = mkpp(xx, D2P); D2X = ppval(D2px,ts);
        
        % curvature
        K = (DX.*D2Y-DY.*D2X)./(DX.^2+DY.^2).^1.5;
        
        % speed
        V = sqrt(DX.^2+DY.^2);
    else strcmp(Lcurve_method,'polynomial')
        % not adapted for curvature!
        % splines are better!
        polyfit_order = length(mehat)-1;
        %polyfit_order = 4;
        pme = polyfit(t,mehat,polyfit_order);
        pme_p = polyderiv(pme);
        pme_pp = polyderiv(pme_p);
        X = polyval(pme,ts);
        mepc = polyval(pme_p,ts);
        pre = polyfit(t,rehat,polyfit_order);
        pre_p = polyderiv(pre);
        pre_pp = polyderiv(pre_p);
        Y = polyval(pre,ts);
        repc = polyval(pre_p,ts);
        K = ( mepc.*polyval(pre_pp,ts) - repc.*polyval(pme_pp,ts) ) ./ ( mepc.^2 + repc.^2 ).^(3/2);
        V = sqrt(mepc.^2 + repc.^2);
    end
    
    % Determine characteristic points
    [~,Dc2] = min(X);
    lim1_minME = limr(Dc2);
    printf('Minimum model error at: %g\n',lim1_minME)
    
    [~,Dc3] = min(Y);
    lim1_minRE = limr(Dc3);
    printf('Minimum regularization error at: %g\n',lim1_minRE)
    
    if Dc2<Dc3
        [~,Dc_K] = max(K(Dc2:Dc3));
        Dc_K = Dc_K + Dc2-1;
    else
        [~,Dc_K] = max(K(Dc3:Dc2));
        Dc_K = Dc_K + Dc3-1;
    end
    lim1_maxcurv = limr(Dc_K);
    printf('Maximum curvature at: %g\n',lim1_maxcurv)
    
    if Dc2<Dc3
        [~,Dc_V] = min(V(Dc2:Dc3));
        Dc_V = Dc_V + Dc2-1;
    else
        [~,Dc_V] = min(V(Dc3:Dc2));
        Dc_V = Dc_V + Dc3-1;
    end
    lim1_minspeed = limr(Dc_V);
    printf('Minimum speed at: %g\n',lim1_minspeed)
    
    %npoints = 2;
    %diff_ratio = 10;
    %pmel = polyfit(ts(Dc2:Dc2+N*(npoints-1)), X(Dc2:Dc2+N*(npoints-1)), 1);
    %Xl = polyval(pmel, ts(Dc2:Dc3));
    %sdXl = std(X(Dc2:Dc2+N*(npoints-1))-Xl(1:1+N*(npoints-1)));
    %Dc_signME = Dc2 -1 + find((X(Dc2:Dc3)-Xl) > diff_ratio*sdXl)(1);
    %lim1_signME = limr(Dc_signME)
    
    % Plot L-curve
    figure(2); plot(mehat,rehat,'rx',X,Y,'b-',X(Dc_K),Y(Dc_K),'gx',X(Dc_V),Y(Dc_V),'mx',X(Dc2),Y(Dc2),'cx')
    xlabel('Log10 ME')
    ylabel('Log10 RE')
    legend('samples','fit','max curvature','min speed','min ME', -1)
    axis('equal')
    title('L-curve')
    
    % Plot normalized curvature and speed
    figure(3);
    Kma = max(K); Kmi = min(K);
    Vma = max(V); Vmi = min(V);
    plot(limr, (K-Kmi)/(Kma-Kmi), 'r', limr, (V-Vmi)/(Vma-Vmi), 'b')
    xlabel('Log10 lim')
    ylabel('Normalized curvature and speed')
    axis('normal')
    legend('Curvature K', 'Speed V')
    title('')
    
    % Save results
    name = sprintf('%s/Lcurve.mat', resultdir);
    save(name, 'lim1_maxcurv', 'lim1_minspeed','lim1_minME', 'lim1_minRE', 'K', 'V', 'me', 're', 'lim1', 'limr')
    end
