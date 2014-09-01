

function [c1] = gen_calls(c0, insensitivity)
    [m,n] = size(c0)
    c1 = c0 | (grand(m,n,"uin", 0, 100) < insensitivity)
endfunction


function [f1,h1, u1,d1,i1] = service_calls(f0,h0,u0,d0,i0)
    f1 = f0
    h1 = h0
    u1 = u0
    d1 = d0
    i1 = i0

    top = size(i0,'c') - 1

    // Are there any calls?
    calls = u0 | d0 | i0
    if ~or(calls) then
        h1 = 0
        return
    end

    // Adjust heading, 
    calls_f = find(calls) - 1
    f_max = max(calls_f)
    f_min = min(calls_f)
    if f0 <= f_min then
        if f0 == f_min then
            h1 = 0
        else
           // no calls below, head up (includes bottom)
            h1 = +1
        end
    elseif f_max <= f0 then
        if f0 == f_max then
            h1 = 0
        else
            // no calls above, head down (includes top)
            h1 = -1
        end
    elseif (h0 == 0) then
        // head toward the closest call 
        if abs(f0 - f_min) < abs(f_max - f0) then
            if  f0 < f_min then
                h1 = +1
            elseif f_min < f0 then
                h1 = 0
            else
                h1 = -1
            end
        else
            if f0 < f_max then
                h1 = +1
            elseif f_max < f0 then
                h1 = 0
            else
                h1 = -1
            end
        end
    else
        // prefer to continue on current heading
        h1 = h0
    end

    // Update the current floor. 
    f1 = f0 + h1

    // Service the floor by clearing its calls.
    fix = f1 + 1
    // mprintf("fix %d %d %d %d %d\n", f0, fix, f_min, f_max, h1)
    u1(fix) = 0
    d1(fix) = 0
    i1(fix) = 0

endfunction

function s = format_vector(v0,v1)
    s = ""
    v2 = v0 & v1
    for ix = 1 : size(v1,'c')
        if v1(ix) then
            s = s + "X"
        elseif v0(ix) & ~v1(ix) then
            s = s + "O"
        else
            s = s + "_"
        end
        s = s + " "
    end
endfunction

function s = format_floor(value, upper)
    s = ""
    for ix = 1 : upper
        if (ix-1) == value then
            s = s + "X"
        else
            s = s + "_"
        end
        s = s + " "
    end
endfunction

function s = format_heading(value)
    if value < 0 then
        s = "<<<<<"
    elseif 0 < value then
        s = ">>>>>"
    else
        s = "====="
    end
endfunction

function [] = display_state(action, t, f,h,u0,u1,d0,d1,i0,i1) 
    mprintf("%s : id(%d)\n", action, t)

    mprintf("%s : %s \n", "Heading", format_heading(h))
    mprintf("%s : %s \n", "Floor  ", format_floor(f, size(u,'c')))
    mprintf("%s : %s \n", "UP     ", format_vector(u0,u1))
    mprintf("%s : %s \n", "DOWN   ", format_vector(d0,d1))
    mprintf("%s : %s \n", "INSIDE ", format_vector(i0,i1))
endfunction

function [] = display_state_vars(action,u,d,i) 
    mprintf("%s\n", action)

    mprintf("%s : %s \n", "UP     ", format_vector(u,u))
    mprintf("%s : %s \n", "DOWN   ", format_vector(d,d))
    mprintf("%s : %s \n", "INSIDE ", format_vector(i,i))
endfunction
