function R = radonmtx(dim, angles, offsetnr)

rows = dim^2;
cols = offsetnr*numel(angles);
R = sparse(rows, cols);
disp('Grab a coffee, this may take a while!')
parfor i=1:cols
    e = zeros(cols, 1);
    e(i) = 1;
    R(:, i) = reshape(iradon(reshape(e, [offsetnr numel(angles)]), ...
        angles, 'linear', 'None', 1, dim), [1 rows]);    
end

R = R';

end