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
If not in-motion elevator will move in the direction where the farthest call is closer.
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


