% This performs a study providing a state to generate new calls.
% 

% A test consisting of a single inside call
u = [0 0 0 0 0 0];
d = [0 0 0 0 0 0];
f = 5;
h = +1;
i = [0 0 0 0 0 1];
display_state('Trial A',0,h,f,u,u2,d,d2,i,i2);
[h,f,u2,d2,i2] = handle_event(h,f,u,d,i);
display_state('Trial A',1,h,f,u,u2,d,d2,i,i2);


% A test consisting of two inside calls where the
% heading is downward and the current floor is
% on the upper call.
f = 5;
h = -1;
i = [1 0 0 0 0 1];
display_state('Trial B',0,h,f,u,u2,d,d2,i,i2);
for ix = 1 : 8
    [h,f,u2,d2,i2] = handle_event(h,f,u,d,i);
    display_state('Trial B',ix,h,f,u,u2,d,d2,i,i2);
    u = u2;
    d = d2;
    i = i2;
end

% A trial making sure that once a heading is
% set it does not change until all calls
% on that heading have been serviced.
f = 1;
h = +1;
i = [0 0 0 0 0 0];
d = [0 1 0 0 0 1];
display_state('Trial C',0,h,f,u,u2,d,d2,i,i2);
for ix = 1 : 8
    [h,f,u2,d2,i2] = handle_event(h,f,u,d,i);
    display_state('Trial C',ix,h,f,u,u2,d,d2,i,i2);
    u = u2;
    d = d2;
    i = i2;
end


% This trial makes calls at a fairly high rate.
display_state('Trial D',0,h,f,u,u2,d,d2,i,i2)
for ix = 1 : 20
    [h,f,u2,d2,i2] = handle_event(h,f,u,d,i);
    % display the post-service matrix
    display_state('Trial D',ix,h,f,u,u2,d,d2,i,i2)
    % make calls (randomly)
    u = gen_calls(u2,10);
    d = gen_calls(d2,10);
    i = gen_calls(i2,10);
end

% A trial where the calls are made infrequently.
% This is a continuation of the previous trial.
for ix = 1 : 20
    [h,f,u2,d2,i2] = handle_event(h,f,u,d,i);
    % display the post-service matrix
    display_state('Trial D+',ix,h,f,u,u2,d,d2,i,i2)
    % make calls (randomly)
    u = gen_calls(u2,2);
    d = gen_calls(d2,2);
    i = gen_calls(i2,2);
end
