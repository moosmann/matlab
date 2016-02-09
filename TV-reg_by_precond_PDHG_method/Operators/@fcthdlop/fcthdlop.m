classdef fcthdlop < linop   
    
    properties(SetAccess = private, GetAccess = private)
        fwfcthdl
        bwfcthdl
    end
    
    methods(Access = public)
        function obj = fcthdlop(domaindim, imagedim, fwfcthdl, bwfcthdl)
            obj.prop.type = {'fcthdlop'};
            if (nargin == 0)
                obj.domaindim = [1, 1];
                obj.imagedim = [1, 1];
                obj.fwfcthdl = @(x) x;
                obj.bwfcthdl = @(x) x;
            elseif (nargin == 4)
                obj.domaindim = domaindim;
                obj.imagedim = imagedim;
                obj.fwfcthdl = fwfcthdl;
                obj.bwfcthdl = bwfcthdl;
            else
                error('Not enough input arguments!')
            end
        end
    end
    
    methods(Access = protected)
        function u = backwardmult(obj, f)
            u = obj.bwfcthdl(f);
        end
        function f = forwardmult(obj, u)
            f = obj.fwfcthdl(u);
        end
    end
    
end

