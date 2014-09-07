
function [h1,f1, u1,d1,i1] = handle_event(h0,f0,u0,d0,i0)
%SERVICE_CALLS Service the elevator calls
%
[h1] = update_heading(h0,f0,u0,d0,i0);

% visit the next floor
f1 = f0 + h1;

% Service the floor by clearing its calls.
[u1,d1,i1] = service_calls(h1,f1,u0,d0,i0);

end

