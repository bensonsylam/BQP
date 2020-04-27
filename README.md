# Binary Quadratic Programming (BQP)
A solver for BQP problem under three settings: unconstrained, linear equality constrainted and linearly inequality constrained problems. 

The general BQP problem is

max s<sup>T</sup>Qs s.t. s∈{-1,1}<sup>n</sup>, s∈Ω,with Ω={s|A<sub>ineq</sub> s≤b<sub>ineq</sub>,A<sub>eq</sub> s=b<sub>eq</sub> }

Q is a symmetric matrix. A_ineq and A_eq are the coefficients of the inequality and equality problems. b_ineq and b_eq are the right-hand-side values of the inequality and equality problems

The Matlab implementation can be found in this hub.

