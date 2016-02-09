classdef linop
    
    properties(SetAccess = protected, GetAccess = protected)
        arg1
        arg2
        domaindim
        flag = 'regular';        
        imagedim
        prop = struct('transp', false);                
    end       
    
    methods(Access = public)
        function obj = linop(arg1, arg2, flag)
            obj.domaindim = [1, 1];
            obj.imagedim = [1, 1];
            if (nargin == 1)
                obj.domaindim = arg1;
                obj.imagedim = arg1;
            elseif (nargin == 3)
                obj.flag = flag;
                switch obj.flag
                    case 'scalmult'
                        obj.arg1 = arg1;
                        obj.arg2 = arg2;
                        obj.domaindim = sizend(arg1, 2);
                        obj.imagedim = sizend(arg1, 1);
                        obj.flag = 'mult';
                    case 'matmult'
                        if (isequal(size(arg1, 2), size(arg2, 1)))
                            obj.arg1 = arg1;
                            obj.arg2 = arg2;
                            obj.domaindim = [size(arg2, 2) 1];
                            obj.imagedim = [size(arg1, 1) 1];
                        else
                            error('Dimensions do not match!')
                        end
                        obj.flag = 'mult';
                    case 'scaladd'
                        obj.arg1 = arg1;
                        obj.arg2 = arg2;
                        obj.domaindim = sizend(arg1, 2);
                        obj.imagedim = sizend(arg1, 1);
                    case 'matadd'
                        if (isequal(size(arg1, 2), size(arg2, 2)) && ...
                                isequal(size(arg1, 1), size(arg2, 1)))
                            obj.arg1 = arg1;
                            obj.arg2 = arg2;
                            obj.domaindim = [size(arg1, 2) 1];
                            obj.imagedim = [size(arg1, 1) 1];
                        else
                            error('Dimensions do not match!')
                        end                    
                    case 'mathorzcat'
                        if isequal(size(arg1, 1), size(arg2, 1))
                            obj.arg1 = arg1;
                            obj.arg2 = arg2;
                            obj.domaindim = [(size(arg1, 2) + size(arg2,...
                                2)) 1];
                            obj.imagedim = [size(arg1, 1) 1];
                        else
                            error('Dimensions do not match!')
                        end
                    case 'matvertcat'
                        if isequal(size(arg1, 2), size(arg2, 2))
                            obj.arg1 = arg1;
                            obj.arg2 = arg2;
                            obj.domaindim = [size(arg1, 2) 1];
                            obj.imagedim = [(size(arg1, 1) + size(arg2, ...
                                1)) 1];
                        else
                            error('Dimensions do not match!')
                        end
                end
            elseif (nargin ~= 0) && (nargin ~= 1) && (nargin ~= 3)
                error('Invalid number of input arguments!')
            end
        end
    end
    
    methods(Access = public, Sealed)
        function obj = ctranspose(obj)
            obj.prop.transp = ~obj.prop.transp;            
            aux = obj.domaindim;
            obj.domaindim = obj.imagedim;
            obj.imagedim = aux;
            if ~strcmp(obj.flag, 'regular')
                obj.arg1 = obj.arg1';
                obj.arg2 = obj.arg2';
            end
        end
        function f = feval(arg1, arg2)
            if isequal(size(arg2), [size(arg1, 2) 1])
                switch arg1.flag         
                    case 'regular'
                        switch checkinputs(arg1, arg2)
                            case 'linear'
                                f = mtimes(arg1, arg2);
                            case 'nonlinear'                                
                                f = reshape(nonlinfcteval(arg1, reshape(...
                                    arg2, sizend(arg1, 2))), [size(arg1,...
                                    1) 1]);
                        end                        
                    case 'mult'
                        if ~arg1.prop.transp
                            if isnumeric(arg1.arg1)
                                f = arg1.arg1*(arg1.arg2.feval(arg2));
                            elseif isnumeric(arg1.arg2)
                                f = arg1.arg1.feval(arg1.arg2*arg2);
                            else
                                f = arg1.arg1.feval(arg1.arg2.feval(arg2));
                            end
                        else
                            if isnumeric(arg1.arg2)
                                f = arg1.arg2*(arg1.arg1.feval(arg2));
                            elseif isnumeric(arg1.arg1)
                                f = arg1.arg2.feval(arg1.arg1*arg2);
                            else
                                f = arg1.arg2.feval(arg1.arg1.feval(arg2));
                            end
                        end    
                    case 'matadd'
                        if isnumeric(arg1.arg1)
                            f = arg1.arg1*arg2 + arg1.arg2.feval(arg2);
                        elseif isnumeric(arg1.arg2)
                            f = arg1.arg1.feval(arg2) + arg1.arg2*arg2;
                        else
                            f = arg1.arg1.feval(arg2) + arg1.arg2.feval(...
                                arg2);
                        end                    
                    case 'scaladd'
                        f = (arg1.arg1.feval(arg2)) + arg1.arg2*(ones(...
                            size(arg1.arg1, 1), 1)*(ones(1, ...
                            size(arg1.arg1, 2))*arg2));
                    case 'mathorzcat'
                        if ~arg1.prop.transp
                            if isnumeric(arg1.arg1)
                                f = arg1.arg1*arg2(1:size(arg1.arg1, 2))...
                                    + arg1.arg2.feval(arg2((size( ...
                                    arg1.arg1, 2) + 1):numel(arg2)));
                            elseif isnumeric(arg1.arg2)
                                f = arg1.arg1.feval(arg2(1:size( ...
                                    arg1.arg1, 2))) + arg1.arg2*arg2((...
                                    size(arg1.arg1, 2) + 1):numel(arg2));
                            else
                                f = arg1.arg1.feval(arg2(1:size(...
                                    arg1.arg1, 2))) + arg1.arg2.feval(...
                                    arg2((size(arg1.arg1, 2) + 1):numel(...
                                    arg2)));
                            end
                        else
                            if isnumeric(arg1.arg1)
                                f = [arg1.arg1*arg2; arg1.arg2.feval(...
                                    arg2)];
                            elseif isnumeric(arg1.arg2)
                                f = [arg1.arg1.feval(arg2); ...
                                    arg1.arg2*arg2];
                            else
                                f = [arg1.arg1.feval(arg2); ...
                                    arg1.arg2.feval(arg2)];
                            end
                        end
                    case 'matvertcat'
                        if ~arg1.prop.transp
                            if isnumeric(arg1.arg1)
                                f = [arg1.arg1*arg2; arg1.arg2.feval(...
                                    arg2)];
                            elseif isnumeric(arg1.arg2)
                                f = [arg1.arg1.feval(arg2); ...
                                    arg1.arg2*arg2];
                            else
                                f = [arg1.arg1.feval(arg2); ...
                                    arg1.arg2.feval(arg2)];
                            end
                        else
                            if isnumeric(arg1.arg1)
                                f = arg1.arg1*arg2(1:size(arg1.arg1, 2))...
                                    + arg1.arg2.feval(arg2((size( ...
                                    arg1.arg1, 2) + 1):numel(arg2)));
                            elseif isnumeric(arg1.arg2)
                                f = arg1.arg1.eval(arg2(1:size( ...
                                    arg1.arg1, 2))) + arg1.arg2*arg2((...
                                    size(arg1.arg1, 2) + 1):numel(arg2));
                            else
                                f = arg1.arg1.feval(arg2(1:size( ...
                                    arg1.arg1, 2))) + arg1.arg2.feval( ...
                                    arg2((size(arg1.arg1, 2) + 1):numel(...
                                    arg2)));
                            end
                        end                    
                end
            else
                error('Dimensions do not match!')
            end
        end
        function prop = getproperties(obj)
            prop = obj.prop;
        end
        function obj = gradient(arg1, arg2)
            if isequal(size(arg2), [size(arg1, 2) 1])
                switch arg1.flag         
                    case 'regular'
                        obj = nonlingradeval(arg1, arg2);                    
                    case 'mult'
                        if ~arg1.prop.transp
                            if isscalar(arg1.arg2) && isnumeric(arg1.arg2)
                                obj = arg1.arg2*(arg1.arg1.gradient(arg2));
                            elseif isnumeric(arg1.arg1)
                                obj = arg1.arg1*(arg1.arg2.gradient(arg2));                               
                            elseif isnumeric(arg1.arg2)                                
                                obj = arg1.arg1.gradient(...
                                    arg1.arg2*arg2)*arg1.arg2;                                
                            else
                                obj = arg1.arg1.gradient(...
                                    arg1.arg2.feval(...
                                    arg2))*arg1.arg2.gradient(arg2);
                            end
                        else
                            if isscalar(arg1.arg2) && isnumeric(arg1.arg2)
                                obj = arg1.arg2*(arg1.arg1.gradient(...
                                    arg2)');
                            elseif isnumeric(arg1.arg1)
                                obj = arg1.arg2.gradient(arg2)'*arg1.arg2';                                
                            elseif isnumeric(arg1.arg2)
                                obj = arg1.arg2'*arg1.arg1.gradient(...
                                    arg1.arg2*arg2);
                            else
                                obj = arg1.arg2.gradient(...
                                    arg2)'*arg1.arg1.gradient(...
                                    arg1.arg2.feval(arg2))';
                            end
                        end    
                    case 'matadd'
                        if isnumeric(arg1.arg1)
                            obj = arg1.arg1 + arg1.arg2.gradient(arg2);
                        elseif isnumeric(arg1.arg2)
                            obj = arg1.arg1.gradient(arg2) + arg1.arg2;
                        else
                            obj = arg1.arg1.gradient(arg2) + ...
                                arg1.arg2.gradient(arg2);
                        end                    
                    case 'scaladd'
                        obj = (arg1.arg1.gradient(arg2)) + arg1.arg2;
                    case 'mathorzcat'
                        if ~arg1.prop.transp
                            if isnumeric(arg1.arg1)
                                obj = [arg1.arg1, arg1.arg2.gradient(...
                                    arg2((size(arg1.arg1, 2) + 1):numel(...
                                    arg2)))];
                            elseif isnumeric(arg1.arg2)
                                obj = [arg1.arg1.gradient(arg2(1:size(...
                                    arg1.arg1, 2))), arg1.arg2];
                            else
                                obj = [arg1.arg1.gradient(arg2(1:size(...
                                    arg1.arg1, 2))), arg1.arg2.gradient(...
                                    arg2((size(arg1.arg1, 2) + 1):numel(...
                                    arg2)))];
                            end
                        else
                            if isnumeric(arg1.arg1)
                                obj = [arg1.arg1; arg1.arg2.gradient(...
                                    arg2)];
                            elseif isnumeric(arg1.arg2)
                                obj = [arg1.arg1.gradient(arg2); ...
                                    arg1.arg2];
                            else
                                obj = [arg1.arg1.gradient(arg2); ...
                                    arg1.arg2.gradient(arg2)];
                            end
                        end
                    case 'matvertcat'
                        if ~arg1.prop.transp
                            if isnumeric(arg1.arg1)
                                obj = [arg1.arg1; arg1.arg2.gradient(...
                                    arg2)];
                            elseif isnumeric(arg1.arg2)
                                obj = [arg1.arg1.gradient(arg2); ...
                                    arg1.arg2];
                            else
                                obj = [arg1.arg1.gradient(arg2); ...
                                    arg1.arg2.gradient(arg2)];
                            end
                        else
                            if isnumeric(arg1.arg1)
                                obj = [arg1.arg1, arg1.arg2.gradient(...
                                    arg2((size(arg1.arg1, 2) + 1):numel(...
                                    arg2)))];
                            elseif isnumeric(arg1.arg2)
                                obj = [arg1.arg1.gradient(arg2(1:size(...
                                    arg1.arg1, 2))), arg1.arg2];
                            else
                                obj = [arg1.arg1.gradient(arg2(1:size(...
                                    arg1.arg1, 2))), arg1.arg2.gradient(...
                                    arg2(size(arg1.arg1, 2) + 1):numel(...
                                    arg2))];
                            end
                        end                    
                end
            else
                error('Dimensions do not match!')
            end            
        end
        function obj = horzcat(arg1, arg2, varargin)
            if (nargin > 2) && (numel(varargin) > 1)
                obj = horzcat(horzcat(arg1, arg2), varargin{1}, ...
                    varargin{2:numel(varargin)});
            elseif (nargin > 2)
                obj = horzcat(horzcat(arg1, arg2), varargin{1});
            else
                switch checkinputs(arg1, arg2)
                    case 'linear'
                        obj = linop(arg1, arg2, 'mathorzcat');
                    case 'nonlinear'
                        obj = nonlinop(arg1, arg2, 'mathorzcat');
                    otherwise
                        error(['Objects of type linop/nonlinop can ' ...
                            'only be concatenated with matrices, or ' ...
                            'objects of the same type!'])                    
                end
            end
        end
        function obj = minus(arg1, arg2)
            if isnumeric(arg1)
                obj = minus(arg2, arg1);
            else
                obj = plus(arg1, -arg2);
            end
        end
        function f = mtimes(arg1, arg2)
            if isnumeric(arg1) && isscalar(arg1)
                f = mtimes(arg2, arg1);
            else
                if isscalar(arg2)
                    switch checkinputs(arg1, arg2)
                        case 'linear'
                            f = linop(arg1, arg2, 'scalmult');
                        case 'nonlinear'
                            f = nonlinop(arg1, arg2, 'scalmult');
                    end
                elseif isnumeric(arg2) && isvector(arg2) && isequal( ...
                        size(arg2, 2), 1)
                    switch arg1.flag
                        case 'regular'
                            if ~arg1.prop.transp
                                f = reshape(forwardmult(arg1, reshape( ...
                                    arg2, sizend(arg1, 2))), [size(arg1,...
                                    1) 1]);
                            else
                                f = reshape(backwardmult(arg1, reshape( ...
                                    arg2, sizend(arg1, 2))), [size(arg1,...
                                    1) 1]);
                            end
                        case 'mult'
                            if ~arg1.prop.transp
                                f = arg1.arg1*(arg1.arg2*arg2);
                            else
                                f = arg1.arg2*(arg1.arg1*arg2);
                            end
                        case 'matadd'
                            f = arg1.arg1*arg2 + arg1.arg2*arg2;
                        case 'scaladd'
                            f = (arg1.arg1*arg2) + arg1.arg2*(ones(...
                                size(arg1.arg1, 1), 1)*(ones(1, ...
                                size(arg1.arg1, 2))*arg2));
                        case 'mathorzcat'
                            if ~arg1.prop.transp
                                f = arg1.arg1*arg2(1:size(arg1.arg1, 2))...
                                    + arg1.arg2*arg2((size(arg1.arg1, 2)...
                                    + 1):numel(arg2));
                            else
                                f = [arg1.arg1*arg2; arg1.arg2*arg2];
                            end
                        case 'matvertcat'
                            if ~arg1.prop.transp
                                f = [arg1.arg1*arg2; arg1.arg2*arg2];
                            else
                                f = arg1.arg1*arg2(1:size(arg1.arg1, 2))...
                                    + arg1.arg2*arg2((size(arg1.arg1, 2)...
                                    + 1):numel(arg2));
                            end
                    end                
                elseif isequal(size(arg1, 2), size(arg2, 1))
                    switch checkinputs(arg1, arg2)
                        case 'linear'
                            f = linop(arg1, arg2, 'matmult');
                        case 'nonlinear'
                            f = nonlinop(arg1, arg2, 'matmult');
                    end
                else
                    error('Dimensions do not match!')
                end
            end
        end
        function val = normest(obj, varargin)
            maxiter = 100;
            tol = 10^-5;
            if numel(varargin) >= 1
                maxiter = varargin{1};
            elseif numel(varargin) == 2
                tol = varargin{2};
            elseif numel(varargin) > 2
                error('Too many input arguments!')
            end
            u = randn([size(obj, 2) 1]);
            u = u/norm(u);
            k = 1;
            sens = Inf;
            while (k <= maxiter) && (sens > tol)
                valold = norm(obj*u)/norm(u);
                u = norm(u)/(norm(obj*u)^2) * (obj'*(obj*u));
                val = norm(obj*u)/norm(u);
                sens = abs(val - valold)/val;
            end            
        end
        function obj = plus(arg1, arg2)
            if isnumeric(arg1)
                obj = plus(arg2, arg1);
            else
                if isscalar(arg2)
                    switch checkinputs(arg1, arg2)
                        case 'linear'
                            obj = linop(arg1, arg2, 'scaladd');
                        case 'nonlinear'
                            obj = nonlinop(arg1, arg2, 'scaladd');                        
                    end
                else
                    switch checkinputs(arg1, arg2)
                        case 'linear'
                            obj = linop(arg1, arg2, 'matadd');
                        case 'nonlinear'
                            obj = nonlinop(arg1, arg2, 'matadd');
                        otherwise
                            error(['Objects of class linop/nonlinop ' ...
                                'can only be added to matrices or ' ...
                                'objects of class linop!']);
                    end
                end
            end
        end
        function dim = size(obj, arg)
            if nargin == 1
                dim = [prod(obj.imagedim), prod(obj.domaindim)];
            elseif nargin > 2
                error('Too many input arguments!')
            elseif (arg == 1)
                dim = prod(obj.imagedim);
            elseif (arg == 2)
                dim = prod(obj.domaindim);
            end
        end
        function dim = sizend(obj, arg)
            if nargin == 1
                dim = {obj.imagedim, obj.domaindim};
            elseif nargin > 2
                error('Too many input arguments!')
            elseif (arg == 1)
                dim = obj.imagedim;
            elseif (arg == 2)
                dim = obj.domaindim;
            end
        end
        function obj = uminus(obj)
            obj = obj*(-1);
        end
        function obj = uplus(obj)
        end
        function obj = vertcat(arg1, arg2, varargin)
            if (nargin > 2) && (numel(varargin) > 1)
                obj = vertcat(vertcat(arg1, arg2), varargin{1}, ...
                    varargin{2:numel(varargin)});
            elseif (nargin > 2)
                obj = vertcat(vertcat(arg1, arg2), varargin{1});
            else
                switch checkinputs(arg1, arg2)
                    case 'linear'
                        obj = linop(arg1, arg2, 'matvertcat');
                    case 'nonlinear'
                        obj = nonlinop(arg1, arg2, 'matvertcat');
                    otherwise
                        error(['Objects of type linop can only be ' ...
                            'concatenated with matrices or objects of ' ...
                            'type linop/nonlinop!'])
                end
            end
        end
    end
    
    methods(Access = protected)
        function u = backwardmult(~, f)
            u = f;
        end
        function aux = checkinputs(arg1, arg2)
            aux1 = (any(strncmp('nonlinop', superclasses(arg1), 8)) ...
                || strncmp('nonlinop', class(arg1), 8));
            aux2 = (any(strncmp('nonlinop', superclasses(arg2), 8)) ...
                || strncmp('nonlinop', class(arg2), 8));            
            if aux1 || aux2
                aux = 'nonlinear';
            else
                aux = 'linear';
            end
        end
        function f = forwardmult(~, u)
            f = u;
        end      
        function f = nonlinfcteval(arg1, u)
            f = forwardmult(arg1, u);
        end
        function obj = nonlingradeval(arg1, ~)
            obj = arg1;
        end
    end
    
end