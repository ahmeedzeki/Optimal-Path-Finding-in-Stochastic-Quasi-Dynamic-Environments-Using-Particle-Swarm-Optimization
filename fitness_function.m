function [val,obs_t] = fitness_function(x,y,x_target,obs)
%COSTFUNCTION Summary of this function goes here
%   Detailed explanation goes here
% x_target = [45 45];
x_obs = obs.x_obs;
y_obs = obs.y_obs;
w_obs = obs.w_obs;
h_obs = obs.h_obs;
blur = 1;

for j = 1:length(x_obs)
    rx1 = obs.rx*rand(1,1);
    rx2 = obs.rx*rand(1,1);
    ry1 = obs.ry*rand(1,1);
    ry2 = obs.ry*rand(1,1);
    obs_t(j,:) = [x_obs(j)-rx1/2 w_obs(j)+rx2/2 y_obs(j)-ry1/2 h_obs(j)+ry2/2];
    in_obs_x(j,:) = x>obs_t(j,1) & x< (obs_t(j,1)+obs_t(j,2));
    in_obs_y(j,:) = y>obs_t(j,3) & y< (obs_t(j,3)+obs_t(j,4));
end
ind = any(in_obs_x.*in_obs_y);
val(~ind) = (x(~ind)-x_target(1)).^2+(y(~ind)-x_target(2)).^2;
val(ind) = 999999999;
        
end
