clear
clc
%%
iter = 1000;
NP = [2 4 6];
Alpha = [0 100 10000];
Npop = 50:50:500;
% create a table
varNames = {'NP','Alpha','npop','average cost','average rf','var cost','var rf'};
varTypes = {'double','double','double','double','double','double','double'};
table_size = length(NP)*length(Alpha)*length(Npop);
result_table = table('Size',[table_size,7],'VariableTypes',varTypes,'VariableNames',varNames);

table_iter = 0;
wb = waitbar(0,'Please wait!');
for np = NP
    for alpha = Alpha
        for npop = Npop
            table_iter = table_iter+1;
            waitbar(table_iter/table_size)
            file_path = ['results_alpha_' num2str(alpha) '_np_' num2str(np) '\\npop_' num2str(npop) '.mat'];

            load(file_path)

            %% Monte Carlo 
            Avg_cost = zeros(1,iter);
            Avg_rf = Avg_cost;
            Vp = [[x0;y0] GlobalBest.Position' [x_des;y_des]];
            for it = 1:iter
                cst = 0;
                rft = 0;
                    for lm = 1:(length(Vp)-1)
                        [obs_t,cs_temp,rf] = line_integral2_with_risk_factor(Vp(:,lm)',Vp(:,lm+1)',obs);
                        cst = cst + cs_temp;
                        rft = rft + sum(rf(:));

                    end
                    Avg_cost(it) = cst;
                    Avg_rf(it) = rft;

            end
            result_table(table_iter,:) = {np,alpha,npop,mean(Avg_cost),var(Avg_cost),mean(Avg_rf),var(Avg_rf)};
        end
    end
end
close (wb)
writetable(result_table,'result_table')
%% plotting
% result_table = readtable('result_table_obs23');
result_table = readtable('result_table_obs13','PreserveVariableNames',true);
style = ['-bx ';'-.rs';'--ok'];

NP = [2 4 6];
y_axis = result_table.("average_cost");
y = reshape(y_axis,10,[]);
t = 0;
Npop = 50:50:500;
for np = NP
    
        t = t+1;
        % extract variable from table
        x_axis = Npop;
        
        semilogy(x_axis ,y(:,3*t-2),style(1,:),x_axis,y(:,3*t-1),style(2,:),x_axis,y(:,3*t),style(3,:));
        xlabel('Population Number')
        ylabel('Log of Average Cost')
        ylim([.10,1E12])
        legend('\alpha=0','\alpha = 100','\alpha =  10000','Location','southwest')
        SetFigure()
        export_fig ('-transparent', '-pdf', ['MonteCarlo/result_ahmed' '_np_' num2str(np)]) ;
    
end