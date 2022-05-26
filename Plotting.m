%------------ Plotting the result ---------------
%clear;
%load('Final_1000_rep_c');
%load('Final_1000_rep_sc');
p1 = Plot_main(oracle,rel_subopt,[1 0 0],1,numsim);
if sc>0
    p2 = Plot_main(oracle_sc,rel_subopt_sc,[0 1 1],1,numsim);
end
grid 'on'
xlabel({'Number of gradient computations'},'Interpreter','latex');
ylabel({'$|\rho(x)-\rho(x^*)|/|\rho(x^*)|$'},'Interpreter','latex');
if sc>0
    legend([p1 p2],{'APDB1','APDB2'});
else
    legend([p1],{'APDB1'});
end

p1 = Plot_main(oracle,rel_infeas,[1 0 0],2,numsim);
if sc>0
    p2 = Plot_main(oracle_sc,rel_subopt_sc,[0 1 1],2,numsim);
end
grid 'on'
xlabel({'Number of gradient computations'},'Interpreter','latex');
ylabel({'$\frac{1}{n}\sum_{i=1}^n \max\{G_i(x),0\}$'},'Interpreter','latex');
if sc>0
    legend([p1 p2],{'APDB1','APDB2'});
else
    legend([p1],{'APDB1'});
end

function y = Plot_main(oracle,rel_subopt,color,fignum,numsim)
    cell_leng_oracle = cellfun(@length,oracle,'uni',false);
    cell_maxiter_oracle = cellfun(@max,oracle,'uni',false);
    x_axis = linspace(1, max(cell2mat(cell_maxiter_oracle)),max(cell2mat(cell_leng_oracle)))';
    leng_oracle_max = 1:max(cell2mat(cell_leng_oracle));
    for sim=1:numsim
        int_aux = interp1(oracle{1,sim}, rel_subopt{1,sim}, x_axis, 'nearest', 'extrap');
        ind_f(sim,1) = find(int_aux(40:length(int_aux))>0,1)+40;
        int_rel_subopt{1,sim} = int_aux;
    end
    ind_f_min = min(ind_f);
    for sim=1:numsim
        int_rel_subopt{1,sim}(ind_f_min:ind_f(sim,1)) = int_rel_subopt{1,sim}(ind_f(sim,1));
    end
    color_line = color;
    figure(fignum);
    y = semilogy(x_axis, mean(cell2mat(int_rel_subopt),2), 'color', ...
        color_line, 'LineWidth', 1.5);
    hold 'on'
end