
function [h1] = update_heading(h0,f0,u0,d0,i0)
%UPDATE_HEADING Service the elevator calls
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
    if f0 == f_min
        % handle the minimum (ground floor)
        h1 = 0;
    elseif f0 == f_max
        % handle the top floor
        h1 = 0;
    elseif f0 < f_min
        % all calls are above
        h1 = +1;
    elseif (f_max < f0)
        % all calls are below
        h1 = -1;
    else
        % head toward the most remote call
        if (f0 - f_min) < (f_max - f0)
            h1 = +1;
        else
            h1 = -1;
        end
    end
else
    % otherwise continue on current heading
    h1 = h0;
end

end

