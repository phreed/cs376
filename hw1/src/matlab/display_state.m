
function [] = display_state(action, t, f,h,u0,u1,d0,d1,i0,i1)
fprintf('%s : id(%d)\n', action, t)
fprintf('%s : %s \n', 'Heading', format_heading(h))
fprintf('%s : %s \n', 'Floor  ', format_floor(f, size(u0,1)))
fprintf('%s : %s \n', 'UP     ', format_vector(u0,u1))
fprintf('%s : %s \n', 'DOWN   ', format_vector(d0,d1))
fprintf('%s : %s \n', 'INSIDE ', format_vector(i0,i1))
end


function [s] = format_vector(v0,v1)
s = '';
v2 = v0 & v1;
[rows,cols] = size(v1);
for ix = 1 : cols
    if v1(ix)
        s = s + 'X';
    elseif v0(ix) && ~v1(ix)
        s = s + 'O';
    else
        s = s + '_';
    end
    s = s + ' ';
end
end


function [s] = format_floor(value, upper)
s = '';
for ix = 1 : upper
    if (ix-1) == value
        s = s + 'X';
    else
        s = s + '_';
    end
    s = s + ' ';
end
end

function [s] = format_heading(value)
if value < 0
    s = '<<<<<';
elseif 0 < value
    s = '>>>>>';
else
    s = '=====';
end
end
