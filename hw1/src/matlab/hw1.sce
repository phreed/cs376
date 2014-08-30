
u = [0 0 0 0 0 0];
d = [0 0 0 0 0 0];
i = [0 0 0 0 0 0];

f = 0;
h = 0; // down : negative, wait : 0; up : positive

for ix = 1 : 20
    mprintf("\n")
    [f,h,u,d,i] = service_calls(f,h,u,d,i)
    // display the post-service matrix
    display_state("SERVICED-CALLS",ix, f,h,u,d,i)
    
    mprintf("\n")
    // make calls (randomly)
    [u,d,i] = make_calls(u,d,i);
    display_state("MEW-CALLS",ix, f,h,u,d,i)
    mprintf("\n")

end
