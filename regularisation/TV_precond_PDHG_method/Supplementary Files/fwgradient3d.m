function grad = fwgradient3d(u)

grad = zeros([size(u) 3]);
grad(:, 1:(end - 1), :, 1) = u(:, 2:end, :) - u(:, 1:(end - 1), :);
grad(1:(end - 1), :, :, 2) = u(2:end, :, :) - u(1:(end - 1), :, :);
grad(:, :, 1:(end - 1), 3) = u(:, :, 2:end) - u(:, :, 1:(end - 1));

end