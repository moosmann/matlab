function u = bwgradient3d(grad)

dim = size(grad);
u = zeros(dim(1:3));

u(:, 1, :) = -grad(:, 1, :, 1);
u(:, 2:(end - 1), :) = grad(:, 1:(end - 2), :, 1) - grad(:, 2:(end - 1),...
    :, 1);
u(:, end, :) = grad(:, end - 1, :, 1);

u(1, :, :) = u(1, :, :) - grad(1, :, :, 2);
u(2:(end - 1), :, :) = u(2:(end - 1), :, :) + grad(1:(end - 2), :, :, 2)...
    - grad(2:(end - 1), :, :, 2);
u(end, :, :) = u(end, :, :) + grad(end - 1, :, :, 2);

u(:, :, 1) = u(:, :, 1) - grad(:, :, 1, 3);
u(:, :, 2:(end - 1)) = u(:, :, 2:(end - 1)) + grad(:, :, 1:(end - 2), 3)...
    - grad(:, :, 2:(end - 1), 3);
u(:, :, end) = u(:, :, end) + grad(:, :, end - 1, 3);

end