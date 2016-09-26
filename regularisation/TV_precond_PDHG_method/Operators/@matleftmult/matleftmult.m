classdef matleftmult < fcthdlop            
    
    methods(Access = public)
        function obj = matleftmult(M, varargin)
            domaindim = size(M);
            if (nargin > 1) && (numel(varargin{1}) == 2)
                domaindim = varargin{1};
                if ~isequal(size(M, 2), domaindim(1))
                    error('Dimensions do not match!')
                end
            end
            imagedim = [size(M, 1), domaindim(2)];
            fwfcthdl = @(u) M*u;
            bwfcthdl = @(f) M'*f;
            obj = obj@fcthdlop(domaindim, imagedim, fwfcthdl, bwfcthdl);
        end
    end
    
end