function CloseToSingular

A=[1 2 3;...
        4 5 6;...
        7 8 9];
b=A*[1;2;3];
b,
rcond(A),
y=A\b % x should be [1;2;3]r