
function [u1,d1,i1] = service_calls(h1,f1,u0,d0,i0)
%SERVICE_CALLS Service the elevator calls
% Service the calls that are on the current floor.
% If there are no calls on the current heading then
% pick up any calls for that heading.

calls = u0 | d0 | i0;
calls_f = find(calls) - 1;
f_max = max(calls_f);
f_min = min(calls_f);

i1 = i0;
u1 = u0;
d1 = d0;

fix = f1 + 1;
i1(fix) = 0;

if h1 < 0
    if f1 <= f_min
        u1(fix) = 0;
    end
    d1(fix) = 0;
elseif 0 < h1
    u1(fix) = 0;
    if f_max <= f1
        d1(fix) = 0;
    end
else
    u1(fix) = 0;
    d1(fix) = 0;
end

end

