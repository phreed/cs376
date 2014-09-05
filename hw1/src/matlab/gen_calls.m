 

function [c1] = gen_calls(c0, insensitivity)
[m,n] = size(c0);
c1 = c0 | (random('Uniform',0,100,m,n) < insensitivity);
end
