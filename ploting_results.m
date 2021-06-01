clear
clc
%%

NP = 6; % the number of intermediate points
alpha = 10000; % risk factor trade-off constant
File_path = ['results_alpha_',num2str(alpha),'_np_',num2str(NP),'\\','results.txt']; % File path
File_path_pdf = ['results_alpha_',num2str(alpha),'_np_',num2str(NP),'\\','results.pdf']; % File path
fid = fopen(File_path,'r'); % open the file
if fid == -1    % check if the file exists
    display('File not exists!');
    return
end
A = fscanf(fid,'%s');
fclose(fid);

n_pop_loc = strfind(A,'n_pop');
best_cost_loc = strfind(A,'BestCost');
risk_factor_loc = strfind(A,'RiskFactor');
% Extracting n_pop values
t = 1;
for i = 1:length(n_pop_loc)
    n_pop(i) = str2double(A(n_pop_loc(i)+6:best_cost_loc(i)-1));
    best_cost(i) = str2double(A(best_cost_loc(i)+9:risk_factor_loc(i)-2));
    
    if i == length(n_pop_loc)
        risk_factor(i) = str2double(A(best_cost_loc(i)+11:end));
    else
        risk_factor(i) = str2double(A(best_cost_loc(i)+11:n_pop_loc(i+1)-1));
    end
end

%% plotting the data
% make sure the n_pop is ordered

[n_pop_o,I] = sort(n_pop);
best_cost_o = best_cost(I);
risk_factor_o = risk_factor(I);

% plotting the results
h = figure
plot(n_pop_o,best_cost_o);
xlabel('Population')
ylabel('Best Cost')
SetFigure()
export_fig ('-transparent', '-pdf', File_path_pdf)
close (h)