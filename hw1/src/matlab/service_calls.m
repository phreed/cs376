
function [u1,d1,i1] = service_calls(h1,f1,u0,d0,i0)
%SERVICE_CALLS Service the elevator calls
%

i1 = i0;
u1 = u0;
d1 = d0;

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

