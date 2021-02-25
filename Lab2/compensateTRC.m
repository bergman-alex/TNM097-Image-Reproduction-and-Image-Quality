function [out] = compensateTRC(r, g, b, in)

x = 0:0.01:1;
out(:, :, 1) = interp1(r, x, in(:, :, 1), 'pchip');
out(:, :, 2) = interp1(g, x, in(:, :, 2), 'pchip');
out(:, :, 3) = interp1(b, x, in(:, :, 3), 'pchip');

end

