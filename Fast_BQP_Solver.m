%%%%%This m-file is an implementation of the work "A Fast Binary Quadratic Programming Solver based on Stochastic Neighborhood Search " (IEEE TPAMI)
%%%%%Authors: Benson Lam (The Hang Seng University of Hong Kong, Hong Kong) and Alan Liew (Griffith University, Australia)
%%%%%This function attempts to solve the following problem
%%%%% max x'*Q*x
%%%%%s.t. A_ineq*x<=b_ineq
%%%%%s.t. A_eq*x=b_eq
%%%%%x is either -1 or 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%[x,obj_val,L]=Fast_BQP_solver(Q,A_eq,b_eq,A_ineq,b_ineq)
%%%%%Input:
%%%%%Q is a symmetric square matrix
%%%%%A_eq, b_eq are the coefficients for equality constraints problems. If not available, set A_eq = [],b_eq = []
%%%%%A_ineq, b_ineq are the coefficients for inequality constraints problems. If not available, set A_ineq = [],b_ineq = []
%%%%%Output: 
%%%%%x is the resultant binary vector; 
%%%%%obj_val: objective function value, i.e. x'*Q*x
%%%%%L: lagrange multipliers
function [x,obj_val,L]=Fast_BQP_solver(Q,A_eq,b_eq,A_ineq,b_ineq)
b_eq = b_eq(:); b_ineq = b_ineq(:);
%%%Parameters Setting
mm = 5; %number of different trials for the time-step in gradient ascent search
n_repeated = 50; %number of times for bootstrapping

%%%Initialization
Q0 = Q; Q(1:1+size(Q,1):end) = 0; Q = Q/sqrt(sum(Q(:).^2));

%%%Find the initial guess
if isempty([b_eq;b_ineq]) == 1 %Unconstrained caes
    norm_row_Q = sum(abs(Q0)); norm_row_A0 = norm_row_Q'; norm_row_Q = (norm_row_Q(:)).^(-1/2); Qp0 = diag(norm_row_Q)*Q0*diag(norm_row_Q);
    [V,D] = eigs(Qp0 + Qp0',1,'LA'); new_x = sign(V); obj_val = new_x'*Q0*new_x; L = 0;
end

if isempty([b_eq;b_ineq]) == 0 %Constrained case
    A = [A_eq A_ineq];
    %%%Determine the const for the Lagrange constraints
    [V,~] = eigs(Q0 + Q0',1,'LA'); sV = sign(V); 
    if length(b_ineq) > 0 %%%inequality case
        d = max(b_ineq' - (sV'*A_ineq),0);
        dt = max(b_ineq' - (-sV'*A_ineq),0);
        if norm(dt) < norm(d)
            d = dt; V = -V; sV = -sV;
        end;
        d_ineq = d(:);
    else
        d_eq = zeros(length(b_eq),1);
    end
    if length(b_eq) == 0  d_eq = []; end
    if length(b_ineq) == 0  d_ineq = []; end
    d = [d_ineq;d_eq];
    b_val=[b_eq;b_ineq]-d;
    const = sum(abs(Q(:)))/(norm(sV'*A-b_val')^2); 
    
    A2 = A;
    Q1 = [[-sum(b_val.^2) (A2*b_val)'];[(A2*b_val) (Q0/const(1)- A2*A2')]]; 
    norm_row_Q = sum(abs(Q1)); norm_row_Q(norm_row_Q==0) = 1;
    norm_row_Q = (norm_row_Q(:)).^(-1/2); Qp1 = diag(norm_row_Q)*Q1*diag(norm_row_Q);
    
    [V,~] = eigs(Qp1 + Qp1',1,'LA');
    if abs(V(1))< 1e-4 V(1) = 1e-4; end; V = V(2:end)*sign(V(1));
    [obj_val,new_x,L] = sub_solver(Q0,V(:,1),A_eq,b_eq,A_ineq,b_ineq); %Find a feasible solution
end

%%%Optimization
[x,L] = deter_stoch_search(Q0,Q,new_x,A_eq,b_eq,A_ineq,b_ineq,L,mm,n_repeated);

%%%Output the result
obj_val = x'*Q0*x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,L_best] = deter_stoch_search(Q0,Q,new_x,A_eq,b_eq,A_ineq,b_ineq,L_best,mm,n_repeated)
all_soln = [new_x'*(Q0*new_x);new_x]; max_ord = 0; max_obj_val = all_soln(1); max_x = all_soln(2:end); L = L_best;

while 1 == 1
    %%Start the computation
    residual = 1; tolerance = 0; iter = 0;
    while residual > tolerance
        %%%Solve by LM Solver
        x = new_x(:);
        %Update: Deterministic search
        [new_x,obj_val0,L0,change] = binary_projection_update(Q0,Q,mm,x,A_eq,b_eq,A_ineq,b_ineq);
        
        %%%Record the best solution
        if max_obj_val < obj_val0 & change == 1 
            max_obj_val = obj_val0; max_x = new_x; L_best = L0;
        end
        
        %%%Find the residual
        residual = mean(abs(x-new_x)); iter = iter + 1; 
    end
    
    %%%%%%Only keep the best solution
    if max_obj_val > max(all_soln(1,:))
        all_soln = [max_obj_val;max_x];
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%Apply Bootstrapping Technique for Re-initialization
    boost_ind = randi([1 size(Q,1)],size(Q,1),1); weig_vec = histc(boost_ind,[1:size(Q,1)]);
    [valt,ordt] = max(all_soln(1,:)); store_x = all_soln(2:end,ordt);
    all_store_x = store_x.*weig_vec;
    
    %%%%%%Perform bootstrapping and apply normalization to the constraint
    u = weig_vec.*(Q*(weig_vec.*store_x)); v = all_store_x; v(u==0) = 0; nor_const = norm(v)/norm(u); u = u*nor_const;
    [uni_u,ttind1,ttind2] = unique(u);
    
    %Random time-step and random perturbation
    alpha_val = rand(1);temp_val = (1-alpha_val)*u + alpha_val*v;
    temp_temp = temp_val(ttind1);
    indi = randi(length(uni_u),1); indi = indi(1);
    tnew_x = sign(temp_val - temp_temp(indi)); tnew_x(u==0)=0;
    tnew_x = tnew_x(:);    
    if length(b_ineq) > 0 
        if (weig_vec.*tnew_x)'*[A_eq A_ineq]*L_best(:) < 0
            tnew_x = - tnew_x;
        end
    end
    new_x = tnew_x(:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%Count the number of times max_obj_val is the best. 
    if all_soln(1,1) < max_obj_val
        max_ord = 0;
    else
        max_ord = max_ord + 1;
    end
    
    %%%%%%If max_obj_val is the best over n_repeated times, break the loop.
    if max_ord > n_repeated
        break;
    end
end
[val,ord] = max(all_soln(1,:)); x = all_soln(2:end,ord);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_x,obj_val0,L0,change] = binary_projection_update(Q0,Q,mm,x,A_eq,b_eq,A_ineq,b_ineq);
u = Q*x; v = x; nor_const = norm(x)/norm(u); u = u*nor_const; obj_val0 = x'*Q0*x; new_x = x; L0 = 0;
change = 0; 
for alpha_val = 0:1/mm:1-1/mm
    temp_val = (1-alpha_val)*u + alpha_val*v;
    if isempty([b_eq;b_ineq]) == 0 %Constrained caes
        [obj_val,tnew_x,L] = sub_solver(Q0,temp_val,A_eq,b_eq,A_ineq,b_ineq);
    else %Unconstrained caes
        tnew_x = sign(temp_val); obj_val = tnew_x'*(Q0*tnew_x); L = 0;
    end
    if obj_val > obj_val0
        new_x = tnew_x; obj_val0 = obj_val; L0 = L/nor_const; change = 1;
    end
end

