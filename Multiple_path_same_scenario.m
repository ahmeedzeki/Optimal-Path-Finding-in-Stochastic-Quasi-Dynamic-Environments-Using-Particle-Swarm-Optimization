clear
clc
close all
% plot multiple paths on the same scenarion for different alphas

NP = 2;
t = 0;
for alpha = [0 100 10000]
    t = t+1;
    path_to_mat = ['mohanad_BestPath/mat_mohanad/'];
    path_to_mat_file = [path_to_mat 'alpha' num2str(alpha) '_npop500_np' num2str(NP) '.mat'];
    load(path_to_mat_file)
    Points(:,:,t) = BestSol.Position;
end

obs = sim_param.obs;
x0 = sim_param.x0;
y0 = sim_param.y0;
x_des = sim_param.x_des;
y_des = sim_param.y_des;
rec_width = 2;

Vp = [[x0;y0] Points(:,:,1)' [x_des;y_des]];
cst = 0;
rft = 0;
for lm = 1:(length(Vp)-1)
    [obs_t,cs_temp,rf] = line_integral2_with_risk_factor(Vp(:,lm)',Vp(:,lm+1)',obs);
    cst = cst + cs_temp;
    rft = rft + sum(rf(:));

end
style = ['-bx ';'-.rs';'--ok'];

clf()
    hold
%     for j = 1:length(particle)
for j = 1:3
        xax = [x0 Points(:,1,j)' x_des];
        yax = [y0 Points(:,2,j)' y_des];
       plot(xax,yax,style(j,:));
end
legend('\alpha = 0','\alpha = 100','\alpha=10000','Location','SouthEast')
xlabel('x')
ylabel('y')
    rectangle('Position',[x_des-rec_width/2 y_des-rec_width/2 rec_width/2 rec_width/2])
    rectangle('Position',[x0-rec_width/2 y0-rec_width/2 rec_width/2 rec_width/2])
    for j = 1:length(obs.x_obs)
        rectangle('Position',[obs_t(j,1) obs_t(j,3) obs_t(j,2) obs_t(j,4)]);
    end
    for j = 1:length(obs.x_obs)
        rectangle('Position',[obs.x_obs(j)-obs.rx/2 obs.y_obs(j)-obs.ry/2 obs.w_obs(j)+obs.rx obs.h_obs(j)+obs.ry ],'EdgeColor','r');
    end
    xlim([0 50])
    ylim([0 50])
    drawnow 
SetFigure()
export_fig('-pdf','-transparent',[path_to_mat '_npop500_np' num2str(NP) 'simple_scenario'])