
function [f1,h1, u1,d1,i1] = service_calls(f0,h0,u0,d0,i0)
%SERVICE_CALLS Service the elevator calls
%
f1 = f0;
h1 = h0;
u1 = u0;
d1 = d0;
i1 = i0;
top = size(i0,1) - 1;
% Are there any calls?
calls = u0 | d0 | i0;
if not(any(calls))
    h1 = 0;
    return
end
% Adjust heading,
calls_f = find(calls) - 1;
f_max = max(calls_f);
f_min = min(calls_f);
if f0 <= f_min 
    if f0 == f_min 
        h1 = 0;
    else
        % no calls below, head up (includes bottom)
        h1 = +1;
    end
elseif f_max <= f0
    if f0 == f_max
        h1 = 0;
    else
        % no calls above, head down (includes top)
        h1 = -1;
    end
elseif (h0 == 0)
    % head toward the closest call
    if abs(f0 - f_min) < abs(f_max - f0)
        if f0 < f_min
            h1 = +1;
        elseif f_min < f0
            h1 = 0;
        else
            h1 = -1;
        end
    else
        if f0 < f_max
            h1 = +1;
        elseif f_max < f0
            h1 = 0;
        else
            h1 = -1;
        end
    end
else
    % prefer to continue on current heading
    h1 = h0;
end
% Update the current floor.
f1 = f0 + h1;
% Service the floor by clearing its calls.
fix = f1 + 1;
% mprintf("fix %d %d %d %d %d\n", f0, fix, f_min, f_max, h1)
u1(fix) = 0;
d1(fix) = 0;
i1(fix) = 0;
end

