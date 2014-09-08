
function [h1] = adjust_heading(h0,f0,u0,d0,i0)
%ADJUST_HEADING Addust the heading if there are 
% no calls on the current heading.  If there are
% calls on the current heading then keep going.
% If there are no calls then wait.
%

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

% assume we continue on the current heading
if 0 < h0  
    if f0 < f_max
        h1 = h0;
    elseif f0 == f_min
        h1 = 0;
    else
        h1 = -1;
    end
elseif h0 < 0
    if f_min < f0
        h1 = h0;
    elseif f_min == f0
        h1 = 0;
    else
        h1 = +1;
    end
else
    if f_min == f_max && f_min == f0
        h1 = 0;
    elseif f0 < f_min
        h1 = +1;
    elseif f_max < f0
        h1 = -1;
    elseif (f0 - f_min) < (f_max - f0)
        h1 = +1;
    else
        h1 = -1;
    end
end

end

