
function [h1,f1, u1,d1,i1] = handle_event(h0,f0,u0,d0,i0)
%SERVICE_CALLS Service the elevator calls
% This implements the main elevator automata

% Service the floor by clearing its calls.
[u1,d1,i1] = service_calls(h0,f0,u0,d0,i0);

% Update the heading based on the state
[h1] = adjust_heading(h0,f0,u0,d0,i0);

% visit the next floor, where next is identified by the heading
[f1] = visit_floor(h1, f0);

end

