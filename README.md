# Binary Quadratic Programming (BQP)
A solver for BQP problem under three settings: unconstrained, linear equality constrainted and linearly inequality constrained problems. 

The general BQP problem is

max s^T Qs s.t. s∈{-1,1}^n, s∈Ω,with Ω={s|A_ineq s≤b_ineq,A_eq s=b_eq }

Q is a symmetric matrix. A_ineq and A_eq are the coefficients of the inequality and equality problems. b_ineq and b_eq are the right-hand-side values of the inequality and equality problems

The Matlab implementation can be found in this hub.

