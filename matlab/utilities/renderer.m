function [r, d] = renderer(verbose)

if nargin < 1
    verbose = 0;
end

d = opengl('DATA');
r = d.Renderer;

if verbose
    disp(d)
end