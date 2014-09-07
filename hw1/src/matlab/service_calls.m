
function [h1,f1, u1,d1,i1] = service_calls(h0,f0,u0,d0,i0)
%SERVICE_CALLS Service the elevator calls
%
[h1] = update_heading(h0,f0,u0,d0,i0);

% visit the next floor
f1 = f0 + h1;

% Service the floor by clearing its calls.
% fprintf("fix %d %d %d %d %d\n", f0, fix, f_min, f_max, h1) 

u1 = u0;
d1 = d0;
i1 = i0;

fix = f1 + 1;
i1(fix) = 0;

if h1 < 0
    d1(fix) = 0;
elseif 0 < h1
    u1(fix) = 0;
else
    d1(fix) = 0;
    u1(fix) = 0;
end

end

