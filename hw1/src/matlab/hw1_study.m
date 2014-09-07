

% A test consisting of a single inside call
u = [0 0 0 0 0 0];
d = [0 0 0 0 0 0];
f = 5;
h = +1;
i = [0 0 0 0 0 1];
[h,f,u2,d2,i2] = service_calls(h,f,u,d,i);
display_state('Trial C',1,h,f,u,u2,d,d2,i,i2);
fprintf('\n')
 
% A test consisting of two inside calls where the 
% heading is downward and the current floor is  
% on the upper call.
f = 5;
h = -1;
i = [1 0 0 0 0 1];
for ix = 1 : 8
[h,f,u2,d2,i2] = service_calls(h,f,u,d,i);
display_state('Trial D',ix,h,f,u,u2,d,d2,i,i2);
u = u2;
d = d2;
i = i2;
fprintf('\n')
end


% This trial makes calls a a fairly high rate.
u = [0 0 0 0 0 0];
d = [0 0 0 0 0 0];
i = [0 0 0 0 0 0];
f = 0;
h = 0; % down : negative, wait : 0; up : positive
for ix = 1 : 20
    fprintf('\n')
    [h,f,u2,d2,i2] = service_calls(h,f,u,d,i);
    % display the post-service matrix
    display_state('Trial A',ix,h,f,u,u2,d,d2,i,i2)
    fprintf('\n')
    % make calls (randomly)
    u = gen_calls(u2,10);
    d = gen_calls(d2,10);
    i = gen_calls(i2,10);
end

% A trial where the calls are made infrequently. 
% This is a continuation of the previous trial.
for ix = 1 : 20
    fprintf('\n')
    [h,f,u2,d2,i2] = service_calls(h,f,u,d,i);
    % display the post-service matrix
    display_state('Trial A+',ix,h,f,u,u2,d,d2,i,i2)
    fprintf('\n')
    % make calls (randomly)
    u = gen_calls(u2,2);
    d = gen_calls(d2,2);
    i = gen_calls(i2,2);
end
