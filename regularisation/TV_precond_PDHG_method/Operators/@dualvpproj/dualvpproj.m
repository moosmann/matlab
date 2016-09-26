classdef dualvpproj < handle
    
    properties(GetAccess = private, SetAccess = private)        
        dim1
        dim2
        lambda
        proxdata
        sigma1
        sigma2
    end
    
    methods(Access = public)
        
        %% Class constructor
        
        function obj = dualvpproj(imsize, projsize)
            obj.dim1 = prod(projsize);
            obj.dim2 = 3*prod(imsize); 
            obj.lambda = 1;
        end
        
        %% Additional size-method
        
        function dim = size(obj, arg)
            if nargin == 1
                dim = [obj.dim1, obj.dim2];
            elseif nargin > 2
                error('Too many input arguments!')
            elseif (arg == 1)
                dim = obj.dim1;
            elseif (arg == 2)
                dim = obj.dim2;
            end
        end
        
        %% Getter & setter methods
        
        function proxdata = getproxdata(obj)
            proxdata = obj.proxdata(:);
        end
        
        function sigma = getproxparam(obj)
            sigma = spdiags([obj.sigma1; obj.sigma2], 0, sum(size(obj ...
                )), sum(size(obj))); 
        end
        
        function lambda = getregularisationparameter(obj)
            lambda = obj.lambda;
        end
        
        function setproxdata(obj, proxdata)
            obj.proxdata = proxdata(:);
        end
        
        function setproxparam(obj, sigma1, sigma2)
            obj.sigma1 = sigma1;
            obj.sigma2 = sigma2;
        end
        
        function setregularisationparameter(obj, lambda)
            obj.lambda = lambda;
        end
        
        %% Proximity/Resolvent operation
        
        function u = prox(obj, f)
            f1 = f(1:(size(obj, 1)));
            f2 = f((size(obj, 1) + 1):end);
            u1 = (f1 - obj.sigma1.*obj.proxdata)./(1 + obj.sigma1/ ...
                obj.lambda);
            aux = sqrt(sum(abs(reshape(f2, [numel(f2)/3 3])).^2, 2));
            aux = [aux; aux; aux];
            u2 = f2./max(1, aux);
            u = [u1; u2];
        end                                
        
    end
    
end