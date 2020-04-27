# Binary Quadratic Programming (BQP)
A solver for BQP problem under three settings: unconstrained, linear equality constrainted and linearly inequality constrained problems. 

The general BQP problem is

max<sub>s∈{-1,1}<sup>n</sup></sub> s<sup>T</sup>Qs s.t. s∈Ω,with Ω={s|A<sub>ineq</sub> s≤b<sub>ineq</sub>,A<sub>eq</sub> s=b<sub>eq</sub> }

Q is a symmetric matrix. A<sub>ineq</sub> and A<sub>eq</sub> are the coefficients of the inequality and equality constraints. b<sub>ineq</sub> and b<sub>eq</sub> are the right-hand-side values of the inequality and equality constraints.

The Matlab implementation can be found in this hub.

