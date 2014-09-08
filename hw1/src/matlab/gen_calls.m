 

function [c1] = gen_calls(c0, insensitivity)
% GEN_CALLS Generates a set of calls
% Higher insensitivity [1..100] causes a greater
% number of calls to be likely.

[m,n] = size(c0);
c1 = c0 | (random('Uniform',0,100,m,n) < insensitivity);
end
