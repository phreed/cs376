# Introduction to hw1

The project can be simulated using one of the clojure libraries.
[from scratch](http://www.learningclojure.com/2014/01/finite-automata.html)

[using a library](https://github.com/cdorrat/reduce-fsm)

For the homework though I will use Matlab and Simulink.
I don't know if this is the best way.

Part 1:
-------

The following describes my approach and control logic.
This being a model of an elevator certain familiar elevator behavior will be ignored.
For the purpose of this simulation the instantaneous transition from one state to another
can be decomposed into two actions, servicing and moving.
Many familiar activities such as door opening, door closing, holding, boarding or exiting can all be lumpted into the servicing action.
The service action a side effect that includes the clearing of any calls for the prior floor.

In the absense of any active calls the elevator will wait in its current location.
If elevator not in-motion it will move in the direction where the farthest call is closer.
If in-motion elevator will continue in its current direction so long as there are calls in that direction.


Part 2:
-------

The cage heading FSM states:
Moving-Up
Moving-Down
Stationary

The floor FSM has six states, one for each floor

Each button FSM (interior or exterior):
Active 
Inactive

This gives 18 FSM's which compose to form the elevator FSM.
The rule for each FSM are as follows:



f* : no active calls => heading rest
f* : active calls on current floor => clear calls on current floor

f0 : active calls => f1 : clear any f1 calls : heading up
f5 : active calls => f4 : clear any f4 calls : heading down

for i in {1,2,3,4}
fi : heading up : active calls above => f(i+1) : clear f(i+1) calls
fi : no active calls above => f(i+1) : clear f(i+1) calls
fi : heading rest : active calls above => f(i+1) : clear f(i+1) calls : heading up
fi : heading rest : active calls below => f(i-1) : clear f(i-1) calls : heading down
fi : heading down : active calls => f(i-1) : clear f(i-1) calls

