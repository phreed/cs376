
u = [0 0 0 0 0 0];
d = [0 0 0 0 0 0];
i = [0 0 0 0 0 0];
f = 0;
h = 0; % down : negative, wait : 0; up : positive
for ix = 1 : 20
    fprintf('\n')
    [f,h,u2,d2,i2] = service_calls(f,h,u,d,i);
    % display the post-service matrix
    display_state('SERVICED-CALLS',ix, f,h,u,u2,d,d2,i,i2)
    fprintf('\n')
    % make calls (randomly)
    u = gen_calls(u2,10);
    d = gen_calls(d2,10);
    i = gen_calls(i2,10);
end
for ix = 1 : 20
    fprintf('\n')
    [f,h,u2,d2,i2] = service_calls(f,h,u,d,i);
    % display the post-service matrix
    display_state('SERVICED-CALLS',ix, f,h,u,u2,d,d2,i,i2)
    fprintf('\n')
    % make calls (randomly)
    u = gen_calls(u2,2);
    d = gen_calls(d2,2);
    i = gen_calls(i2,2);
end
u = [0 0 0 0 0 0];
d = [0 0 0 0 0 0];
f = 5;
h = +1;
i = [0 0 0 0 0 1];
[f,h,u2,d2,i2] = service_calls(f,h,u,d,i);
display_state('SERVICED-CALLS',1, f,h,u,u2,d,d2,i,i2);
fprintf('\n')
f = 5;
h = -1;
i = [1 0 0 0 0 1];
[f,h,u2,d2,i2] = service_calls(f,h,u,d,i);
display_state('SERVICED-CALLS',1, f,h,u,u2,d,d2,i,i2);
fprintf('\n')