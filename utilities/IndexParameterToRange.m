function x = IndexParameterToRange(x, N)
% Convert paramterized index/indices to a range that can be used for array
% indexing. x can be absolute index values or relative ones. The
% dimension of x can be 1D or 2D. If 1D the start values is x or 1-x and the end
% values is 1 - x or N - x. If x contains relative values integer values
% are computed as x -> round( (N-1) * x + 1 ). N is mandatory for relative
% indexing.
%
% ARGUMENTS
% x: 1D or 2D, float in [-1,1] or integer in [-N,N]. Relative or absolute
% index/indices. if scalar, upper limit is 1-x or N-x, respectively.
% N: scalar, integer. Maximum number of indices. Required for relative
% indexing.
%
% Written by Julian Moosmann, 2016-12. Last version: 2017-02-20
%
% x = IndexParameterToRange(x, N)

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(x,N) round((N - 1) * x + 1 );

numel_x = numel( x );

% Check if x has more than 2 elements
if numel_x > 2
    return
end

% Check if x is empty
if isempty( x )
    x = 1:N;
    return
end

% Check if x has only 1 element
if numel_x == 1
    if x <= -1
        x = [ N + x, N ];
    elseif x > -1 && x < 0
        x = [ 1 + x, 1];
    elseif x >= 0 && x <= 1
        x = [ x, 1 - x ];
    elseif x > 1
        x = [ x, N - x ];
    end
end

% Negative integer indices
x(x<=-1) = x(x<=-1) + N;

% Negative float
m = x > -1 & x < 0;
x(m) = 1 + x(m);

% Convert to range if in [-1 1]

% Check if x(1) is in [0,1)
if x(1) < 1 && x(1) >= 0
    x(1) = f(x(1), N);
end
% Check if x(2) is in (0,1]
if x(2) <= 1 && x(2) > 0
    x(2) = f(x(2), N);
end
x = x(1):x(2);

