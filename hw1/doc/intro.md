# Introduction to hw1

[great documentation](http://jacobian.org/writing/great-documentation/what-to-write/)

The project can be simulated using one of the clojure libraries.
[from scratch](http://www.learningclojure.com/2014/01/finite-automata.html)

[using a library](https://github.com/cdorrat/reduce-fsm)

For the homework though I will use Matlab and Simulink.
I don't know if this is the best way.

Part 1:
-------

In the absense of any active calls the elevator will wait in its current location.
If elevator not in-motion it will move in the direction where the farthest call is closer.
If in-motion elevator will continue in its current direction so long as there are calls in that direction.


Part 2:
-------

The FSM has the following nodes representing the cage heading:
Moving-Up
Moving-Down
Stationary

The FSM has the following state variables:
Floor : f
ExtUp : u x 5
ExtDown : d x 5
IntGoal : i x 6

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

