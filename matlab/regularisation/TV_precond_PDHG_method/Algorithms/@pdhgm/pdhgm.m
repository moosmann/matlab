classdef pdhgm < handle    
    
    properties(GetAccess = private, SetAccess = private)               
        fig
        Fstar
        G
        K
        k
        maxiter
        ploton = false;        
        xprev
        var
    end
    
    methods(Access = public)
        
        %% Class constructor
        
        function obj = pdhgm(K, Fstar, G)
            obj.K = K;
            obj.Fstar = Fstar;
            obj.G = G;
            obj.maxiter = 500;
            obj.k = 1;
        end
        
        %% Additional methods
        
        function disableplot(obj)
            obj.ploton = false;
        end
        
        function enableplot(obj)
            obj.ploton = true;
        end
        
        function initialise(obj)
            obj.k = 1;
            obj.var.x = zeros([size(obj.K, 2) 1]);
            obj.var.y = zeros([size(obj.K, 1) 1]);
        end
        
        function resetcounter(obj)
            obj.k = 1;
        end
        
        %% Setter & getter methods
        
        function maxiter = getmaxiter(obj)
            maxiter = obj.maxiter;
        end
        
        function var = getvariables(obj)
            var = obj.var;
        end       
        
        function setmaxiter(obj, maxiter)
            obj.maxiter = maxiter;
        end       
        
        function setvariables(obj, var)
            obj.var = var;
        end
        
        %% Main routine
        
        function solve(obj)            
            
            if isempty(obj.var)
                obj.initialise;
            end
            
            if obj.ploton
                obj.fig = figure();              
            end
            
            while obj.k <= obj.maxiter
                
                obj.xprev = obj.var.x;
                
                obj.var.x = obj.G.prox(obj.var.x - obj.G.getproxparam * ...
                    (obj.K'*obj.var.y));
                
                obj.var.y = obj.Fstar.prox(obj.var.y + ...
                    obj.Fstar.getproxparam * (obj.K*( 2*obj.var.x - ...
                    obj.xprev )));
                
                if obj.ploton
                    obj.plot
                else
                    disp(['Iteration nr. ' num2str(obj.k) '\' ...
                        num2str(obj.maxiter)])
                end
                
                obj.k = obj.k + 1;
            end
            
            if obj.ploton
                close(obj.fig)
            else
                disp('Computation complete!')
            end
            
        end
                
    end
    
    methods(Access = private)        
        
        %% Method for visualisation
                
        function plot(obj)
            obj.fig;
            aux = reshape(obj.var.x, sizend(obj.G));
            imagesc(aux(:, :, 13))%, [0 1])
            axis image
            colormap(gray(512))
            colorbar
            title(['Iteration nr. ' num2str(obj.k)])
            drawnow
        end
        
    end
    
end