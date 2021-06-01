clear
clc
close all
% plot convergence curves for different alphas same np
MaxIt = 250;
NP = 6;
t = 0;
Best_Cost_all = zeros(MaxIt,3);
scenario = 'hard';


%% Extracting Results
for alpha = [0 100 10000]
    t = t+1;
    path_to_results = ['results_alpha_',num2str(alpha),'_np_',num2str(NP),'\\','results_convergence_',scenario,'.txt'];
    fid = fopen(path_to_results,'r'); % open the file
    if fid == -1    % check if the file exists
        disp('File not exists!');
        return
    end
    A = fscanf(fid,'%s');
    fclose(fid);
    NP_loc = strfind(A,'Np=');
    best_cost_loc = strfind(A,'BestCost=');
    for ii = 1:length(NP_loc)
        if ii == length(NP_loc)
            Best_Cost_all(ii,t) = str2double(A(best_cost_loc(ii)+9:end));
        else
            Best_Cost_all(ii,t) = str2double(A(best_cost_loc(ii)+9:NP_loc(ii+1)-1));
        end
    end

end
%% plotting and exporting figures
File_path_pdf = ['convergence_results_np_',num2str(NP),scenario,'.pdf']; % File path
h = figure
hold on
style = ['-b  ';'-.r ';'--k '];
for j = 1:3
       semilogy(1:MaxIt,Best_Cost_all(:,j),style(j,:));
end
% semilogy(1:MaxIt,Best_Cost_all,style);
xlabel('Iteration')
ylabel('log of Best Cost')
% ylim([0,120000]);
legend('\alpha=0','\alpha=100','\alpha=10000')
SetFigure()
export_fig ('-transparent', '-pdf', File_path_pdf)
close (h)
