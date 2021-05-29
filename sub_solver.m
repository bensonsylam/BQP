%%%%%A function to solve the dual problem
%%%%%This is not the main function
%%%%%This m-file is a function for the "Fast_BQP_Solver.m"
function [obj_val,new_xp,L] = sub_solver(Q,q,A_eq,b_eq,A_ineq,b_ineq)

%%%%When matrices for both inequality and equality constraints are
%%%%non-empty
A = [A_eq A_ineq]; b = [b_eq b_ineq(:)]; lgth_ineq = length(b_ineq); lgth_eq = length(b_eq);
%%%%the lagrange multiplier: L = ru, r is a scalar, u is an unit vector
residual = 1; tolerance = 1e-4; 
tL = zeros(length(b),1);
iter = 0; dL = 1;
sum_A = sum(A);
zero_bd = 0; non_converge = 0;
while residual > tolerance
    iter = iter + 1; tL0 = tL;
    s = sign(q+A*tL); dL = (b - sum_A') + 2*sum(A(s==-1,:))' ;
    if norm(dL)~=0 dLt = dL/norm(dL); [~,tr] = dual_blp( Q,q+A*tL,A*dLt,dLt'*b(:) ); tr = max(tr,0); if lgth_ineq> 0  end; dL = dLt*tr; end;
    if lgth_ineq > 0 
        dL(end-lgth_ineq+1:end) = min(dL(end-lgth_ineq+1:end),zero_bd);
    end
    tL = tL + dL; 
    
    %Find residual
    residual = sum(abs(tL-tL0)./(1+abs(tL0)));
    s = sign(q+A*tL);
    if iter > 100 residual = 0; non_converge = 1; end;
    if lgth_ineq > 0 
        if abs(tr)<tolerance & mean2(max(s'*A-b',0)./sum(A~=0) ) > 0.1 residual = 0; non_converge = 1; end;
    end
    if lgth_eq > 0 
        if abs(s'*A-b')./sum(A~=0) < 0.1
            residual = 0; non_converge = 0;
        end
    end
end

L = tL; new_xp = sign(q+A*L);
if non_converge ==0
    obj_val = new_xp'*Q*new_xp;
else
    obj_val = -inf; new_xp = ones(length(new_xp),1);
end

function [new_x,L] = dual_blp(Q,q,a,b)
ind = a~=0;
if sum(ind)~=0
    [s_val,s_ord] = sort(q(ind)./a(ind)); sorted_a = abs(a(s_ord)); ratio_a = sum(sorted_a) - 2*cumsum(sorted_a);
    [vv,oo] = min( abs(ratio_a - b) ); L = eps*sign(s_val(oo))-s_val(oo); new_x = sign(q+L*a);
else
    L = 0; new_x = sign(q+L*a);
end