function A = fun_times( A, B )
% In-place multiplication of A with B including implicit expansion.

A = B .* A;