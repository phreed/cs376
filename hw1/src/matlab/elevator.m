

function [c1] = gen_calls(c0, insensitivity)
[m,n] = size(c0)
c1 = c0 | (grand(m,n,"uin", 0, 100) < insensitivity)
endfunction
