classdef nonnegproj < handle
    
    properties(GetAccess = private, SetAccess = private)        
        dim
        imsize
        tau
    end
    
    methods(Access = public)
        
        %% Class constructor
        
        function obj = nonnegproj(imsize)
            obj.dim = prod(imsize);
            obj.imsize = imsize;
        end
        
        %% Additional methods
        
        function dim = size(obj)
            dim = obj.dim;
        end
        
        function imsize = sizend(obj)
            imsize = obj.imsize;
        end
        
        %% Getter & setter methods
        
        function tau = getproxparam(obj)
            tau = spdiags(obj.tau, 0, size(obj), size(obj));
        end
        
        function setproxparam(obj, tau)
            obj.tau = tau;
        end
        
        %% Proximity/Resolvent operation
        
        function u = prox(~, f)
            u = max(f, 0);
        end                                
        
    end
    
end