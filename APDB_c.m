function [rel_infeas_err,rel_subopt_err,time_period,iter_epoch,oracle]=...
    APDB_c(input,optval,x0,y0,stop_criteria,max_iter)
%************************************************************
% IMPORTANT: optimal solution is used to measure the accuracy of the
% solution found
%
% Written by Erfan Yazdandoost Hamedani, created on 14 November 2018.
%
% The algorithm is specified to solve Quadratic Constrained 
% Quadratic Programming (QCQP) 
%
% The step-sizes are constant according to 
% Theorem 2.1 part (I) in the paper https://arxiv.org/pdf/1803.01401.pdf
%************************************************************
% min_{x} 0.5*x'*A*x+b'*x
% s.t.    0.5*x'*Q_i*x+d_i'*x-c_i<=0, for i=1:m,
%         -10<=x<=10
%************************************************************
    %---------- Unfolding Input ----------%
    A = input{1,1};
    Q = input{2,1};
    b = input{3,1};
    d = input{4,1};
    c = input{5,1};
    sc = 0;
    [n,~] = size(d);
    optval_rel = max(optval,1);
    %-------------------------------------%
    disp('**********************************************************')
    disp('APDB Convex')
    disp('**********************************************************') 
    epoch=1;
    epoch_counter = 0;
    %------------ step-size Parmeters------%
    eta = 0.7;
    tau = 1e-1;
    gamma = 1;
    sigma_old = gamma*tau;
    tau_old = tau;
    mu = sc;
    %------ Initialization ----------------%
    rel_subopt_err = [];
    rel_infeas_err = [];
    y = y0;
    x = x0;
    iter = 0;
    orc = 0;
    inner = zeros(max_iter,1);
    for j=1:n
       G_y(j,1) = 0.5*x'*Q{j,1}*x+d(j,:)*x-c(j); 
    end
    tic;
    %---------- Main Algorithm ------------%    
    while iter<max_iter
        iter = iter+1;
        G_y_old = G_y;
        for j=1:n
           G_y(j,1) = 0.5*x'*Q{j,1}*x+d(j,:)*x-c(j); 
        end
        orc = orc+1;
        while true
            sigma = gamma*tau;
            theta = sigma_old/sigma;
            ytild = max(y+sigma*((1+theta)*G_y-theta*G_y_old),0);
            G_x = A*x+b;
            for j=1:n
                G_x = G_x+(Q{j,1}*x+d(j,:)')*ytild(j);
            end
            orc = orc+1;
            xtild = max(min(x-tau*G_x,10),-10);
            xdiff = xtild-x;
            E = 0;
            aux = 0;
            for j=1:n
                E = E+((xdiff)'*Q{j,1}*(xdiff))*ytild(j);
                aux = aux+norm(0.5*xtild'*Q{j,1}*xtild-0.5*x'*Q{j,1}*x+d(j,:)*(xdiff))^2;
            end
            E = E+xdiff'*A*xdiff-norm(xdiff)^2/(2*tau)+(sigma)*aux/2;
            if E<=0
                break;
            else
                inner(iter)=inner(iter)+1;
                tau = tau*eta;
            end
        end
        gamma_old = gamma;
        gamma = gamma*(1+mu*tau);
        tau_new = tau*sqrt(gamma_old/gamma *(1+tau/tau_old));
        sigma_old = sigma;
        tau_old = tau;
        tau = tau_new;
        
        y = ytild;
        x = xtild;
        
        if mod(iter,epoch) == 0
            epoch_counter = epoch_counter+1;
            time_period(epoch_counter,1) = toc;
            oracle(epoch_counter,1) = orc;
            subopt = 0.5*x'*A'*x+b'*x;
            infeas = 0;
            for l=1:n
                infeas = infeas+pos(0.5*x'*Q{l,:}*x+d(l,:)*x-c(l));
            end
            rel_subopt_err(epoch_counter,1) = abs(subopt-optval)/abs(optval_rel);
            rel_infeas_err(epoch_counter,1) = infeas/n;
            iter_epoch(epoch_counter,1) = iter;
            if max(abs(subopt-optval)/abs(optval_rel), infeas/n)<stop_criteria
                fprintf(...
                    'Iteration    Time    Rel. Infeas error   Rel. Subopt error\n');
                fprintf('%d    %9.4f       %9.1e         %9.1e\n',iter,time_period(epoch_counter),...
                    rel_infeas_err(epoch_counter),rel_subopt_err(epoch_counter));
                return;
            end
        end
    end
    fprintf(...
        'Iteration    Time    Rel. Infeas error   Rel. Subopt error\n');
    fprintf('%d    %9.4f       %9.1e         %9.1e\n',iter,time_period(epoch_counter),...
        rel_infeas_err(epoch_counter),rel_subopt_err(epoch_counter));
end