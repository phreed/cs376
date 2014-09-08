
function [] = display_state(action, t, h,f,u0,u1,d0,d1,i0,i1)
% DISPLAY_STATE shows the current state of the elevator 
% if the step, t, is zero then this is the initial state.

fprintf('%s : id(%d)\n', action, t)
fprintf('%s : %s \n', 'Heading     ', format_heading(h))
fprintf('%s : %s \n', 'Floor       ', '0 1 2 3 4 5')
fprintf('%s : %s \n', 'Cage Loc    ', format_floor(f, size(u0,2)))
if t < 1
    fprintf('%s : %s \n', 'InCall      ', format_vector(i0,i0))
    fprintf('%s : %s \n', 'ExCall UP   ', format_vector(u0,u0))
    fprintf('%s : %s \n', 'ExCall DOWN ', format_vector(d0,d0))
else
    fprintf('%s : %s \n', 'InCall      ', format_vector(i0,i1))
    fprintf('%s : %s \n', 'ExCall UP   ', format_vector(u0,u1))
    fprintf('%s : %s \n', 'ExCall DOWN ', format_vector(d0,d1))
end
fprintf('\n')
end


function [s] = format_vector(v0,v1)
s = '';
cols = size(v1,2);
for ix = 1 : cols
    if v1(ix)
        s = strcat(s,'X');
    elseif v0(ix) && ~v1(ix)
        s = strcat(s, 'O');
    else
        s = strcat(s,'_');
    end
    s = strcat(s,'.');
end
end


function [s] = format_floor(value, upper)
s = '';
for ix = 1 : upper
    if (ix-1) == value
        s = strcat(s,'X');
    else
        s = strcat(s,'_');
    end
    s = strcat(s,'.');
end
end

function [s] = format_heading(value)
if value < 0
    s = '< < < < < <';
elseif 0 < value
    s = '> > > > > >';
else
    s = '= = = = = =';
end
end
