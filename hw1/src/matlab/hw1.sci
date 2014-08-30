

function [c1] = gen_calls(c0, insensitivity)
    [m,n] = size(c0)
    c1 = c0 | (grand(m,n,"uin", 0, 100) < insensitivity)
endfunction

function [u1,d1,i1] = make_calls(u0,d0,i0) 
    u1 = gen_calls(u0,10)
    d1 = gen_calls(d0,10);
    i1 = gen_calls(i0,10);
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
    if f0 < 1 then
        // bounce at the bottom.
        h1 = +1
    elseif top <= f0
        // bounce at the top
        h1 = -1
    elseif ~(h0 == 0)       
        // prefer the current heading
        h1 = h0
    else  
        // If waiting head toward the farthest call.
        truthy_ix = find(calls) - 1
        if (f0 - min(truthy_ix)) > (max(truthy_ix) - f0) then
            h1 = -1
        else
            h1 = +1
        end 
    end
   
    // Update the current floor. 
    f1 = f0 + h1
    fix = f1+1
    
    // Service the floor by clearing its calls.
    u1(fix) = 0
    d1(fix) = 0
    i1(fix) = 0
  
endfunction

function s = format_vector(vec)
    s = ""
    for ix = 1 : size(vec,'c')
        if vec(ix) then
            s = s + "X"
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

function [] = display_state(action, t, f,h,u,d,i) 
    mprintf("%s : %d\n", action, t)
    
    mprintf("%s : %s \n", "Heading", format_heading(h))
    mprintf("%s : %s \n", "Floor  ", format_floor(f, size(u,'c')))
    mprintf("%s : %s \n", "UP     ", format_vector(u))
    mprintf("%s : %s \n", "DOWN   ", format_vector(d))
    mprintf("%s : %s \n", "INSIDE ", format_vector(i))
endfunction
