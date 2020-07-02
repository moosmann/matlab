function com = CenterOfMass( arr )
% Returns center of mass of an 1-, 2-, or 3-D array.

arr_mean = mean(arr(:));
switch ndims( squeeze( arr ) )
    case {1,2}
    if size( squeeze( arr ) , 1) == 1 || size( squeeze( arr ) , 2) == 1
        x = 1 : max( size( arr ) ); % Columns.
        X = meshgrid(x, 1);
        com = mean(arr(:) .* X(:)) / arr_mean;
        
    else
        x = 1 : size(arr, 2); % Columns.
        y = 1 : size(arr, 1); % Rows.
        [X, Y] = meshgrid(x, y);
        com(2) = mean(arr(:) .* Y(:)) / arr_mean;
        com(1) = mean(arr(:) .* X(:)) / arr_mean;
    end
    case 3
        x = 1 : size(arr, 2); % Columns.
        y = 1 : size(arr, 1); % Rows.
        z = 1 : size(arr, 3); % Rows.
        [X, Y, Z] = meshgrid(x, y, z);
        com(3) = mean(arr(:) .* Z(:)) / arr_mean;
        com(2) = mean(arr(:) .* Y(:)) / arr_mean;
        com(1) = mean(arr(:) .* X(:)) / arr_mean;
end


