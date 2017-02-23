function x = welp()
    x3d = randn(5,6,4)
    new = reshape(permute(x3d,[2 1 3]), 1, [])
end

function [y,z] = fn(n)
    y = 120*n;
    z = 123;
end