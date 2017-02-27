function x = welp()
    x2d = randn(5,6)
    new = reshape(x2d', 1, [])
    transpose(reshape(new, 6, 5))
end
% 
% function [y,z] = fn(n)
%     y = 120*n;
%     z = 123;
% end