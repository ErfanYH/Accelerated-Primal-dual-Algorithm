% Written by Erfan Yazdandoost Hamedani, created on 14 November 2018.
% In this example, APDB algorithm is specified to solve Quadratic Constrained 
% Quadratic Programming (QCQP) 
%
% See section 5.1. of the corresponding paper https://arxiv.org/pdf/1803.01401.pdf
% for more details.
%************************************************************
% min_{x} 0.5*x'*A*x+b'*x
% s.t.    0.5*x'*Q_i*x+d_i'*x-c_i<=0, for i=1:m,
%         -10<=x<=10
%************************************************************
clear;
clc;
seed = 123;
rng(seed,'twister');

numsim = 10;
n = 1e1;
m = 1e3;
x0 = rand(m,1);
y0 = zeros(n,1);
for sim=1:numsim
    S = orth(randn(m,m));
    D = (rand(m-1,1)*100);
    % controlling strong convexity and Lip constants
    A = S'*diag([D;1e-10])*S; % for merely convex scenario
    %A = S'*diag([D;1])*S; % uncomment for strongly convex scenario
    b = randn(m,1);
    sc = min(eig(A));
    if sc<1e-8
        sc = 0;
    end
    d = randn(n,m);
    Q = cell(n,1);
    for j=1:n
        S = orth(randn(m,m));
        D = (rand(m-1,1)*100);
        Q{j,1} = S'*diag([D;1e-10])*S;
    end
    c = rand(n,1);
    %------------- CVX ------------------%
    disp('starting cvx...');
    tic;
    cvx_begin quiet
        cvx_precision high
        variable xstar(m,1)
        dual variable ystar{n}
        constraint = cvx(zeros(n,1));
        for l=1:n
            constraint(l) = 0.5*xstar'*Q{l,:}*xstar+d(l,:)*xstar-c(l);
        end
        minimize(0.5*xstar'*A*xstar+b'*xstar);
        subject to
            for l=1:n
                constraint(l) <= 0 : ystar{l};
            end
            xstar <= 10*ones(m,1);
            xstar >= -10*ones(m,1);
    cvx_end
    time_mosek = toc;
    display(time_mosek);
    optval = cvx_optval;
    %---------------------------------------------%
    input = {A;Q;b;d;c;sc};
    %----------- simulations ---------------------%
    max_iter = 1e5;
    epsilon = 1e-8;

    [rel_infeas{1,sim},rel_subopt{1,sim},time_period{1,sim}...
       ,iter_epoch{1,sim},oracle{1,sim}]=APDB_c(input,optval,x0,y0,epsilon,max_iter);
   
    if sc>0
        [rel_infeas_sc{1,sim},rel_subopt_sc{1,sim},time_period_sc{1,sim}...
             ,iter_epoch_sc{1,sim},oracle_sc{1,sim}]=APDB_sc(input,optval,x0,y0,epsilon,max_iter);
    end
   
end

% Plotting the average of simulations
Plotting;