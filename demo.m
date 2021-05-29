%%%Random Matrix Experiments
%%%Unconstrained problem
clear; load('unconstrained_problem','Q')
[x,obj_val,L]=Fast_BQP_Solver(Q,[],[],[],[]);
display('Random matrix experiement (uniform + normal matrix), n = 1000, section 4.1 of the paper. result shown in appendix IA(ii)')
display('Unconstrained problem')
display(['Objective function value:' num2str(obj_val)])
display('.......................')
%%%Equality constrained problem
clear; load('equality_constrained_problem','Q','A','b')
[x,obj_val,L]=Fast_BQP_Solver(Q,A,b,[],[]);
display('Random matrix experiement (uniform + normal matrix), n = 1000, section 4.1 of the paper. result shown in appendix IB(ii)')
display('Equality constrained problem')
display(['Objective function value:' num2str(obj_val)])
display('.......................')
%%%Inequality constrained problem
clear; load('inequality_constrained_problem','Q','A','b')
[x,obj_val,L]=Fast_BQP_Solver(Q,[],[],A,b);
display('Random matrix experiement (uniform + normal matrix), n = 1000, section 4.1 of the paper. result shown in appendix IC(ii)')
display('Inequality constrained problem')
display(['Objective function value:' num2str(obj_val)])
display('.......................')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image restoration problems:unconstrained problems
% Image restoration problem:square image; noise parameter: 0.1; miu: 0.1; image size: 39x51
load('square data','Q','orig_s');
[x,obj_val,L]=Fast_BQP_Solver(Q,[],[],[],[]); %Apply the proposed BQP solver
figure; imshow(orig_s); title('orig image'); %Show the orignal image
restored_s = reshape(x(1)*x(2:end)>0,39,51); figure; imshow(restored_s); title('restored image'); %Show the restored image using the proposed BQP solver. 
display('Binary image (small) restoration problem')
display('Square image; Noise parameter: 0.1; Miu: 0.1; Image size: 39x51')
display(['Objective function value:' num2str(obj_val)])
display('.......................')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Image restoration problem: rice image; noise parameter: 0.5; miu: 0.5; image size: 256x256
load('big rice data','Q','orig_s');
[x,obj_val,L]=Fast_BQP_Solver(Q,[],[],[],[]); %Apply the proposed BQP solver
figure; imshow(orig_s); title('orig image'); %Show the orignal image
restored_s = reshape(x(1)*x(2:end)>0,256,256); figure; imshow(restored_s); title('restored image'); %Show the restored image using the proposed BQP solver. 
display('Binary image (Big) restoration problem')
display('Rice image; Noise parameter: 0.5; Miu: 0.5; Image size: 256x256')
display(['Objective function value:' num2str(obj_val)])
display('.......................')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image partition problem: equality constrained problems, the constraint is x'1=0
% baby image; sig = 8*255^2; image size: 53x68
load('baby data (ineq partition)','Q','X','A','b');
[x,obj_val,L]=Fast_BQP_Solver(Q,[],[],A,b); %Apply the proposed BQP solver
figure; imshow(X); title('orig image'); %Show the orignal image
partitioned_s = reshape(x>0,53,68); figure; imshow(partitioned_s); title('partition result'); %Show the partitioned image using the proposed BQP solver. 
display('Image partition problem')
display('Baby image; Sig = 8*255^2; Image size: 53x68')
display(['Objective function value:' num2str(obj_val)])
display('.......................')