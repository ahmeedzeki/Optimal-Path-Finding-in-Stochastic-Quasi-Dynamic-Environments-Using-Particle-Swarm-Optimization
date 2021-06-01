function [obs_t,dis,risk] = line_integral2_with_risk_factor(xy1,xy2,obs)
%LINE_INTEGRAL2 Summary of this function goes here
%   Detailed explanation goes here
n_points = 250;
t = linspace(0,1,n_points)';
xyseg = (1-t)*xy1 + t*xy2;
% [xmesh,ymesh] = meshgrid(1:25,1:25);
[zseg,obs_t] = fitness_function(xyseg(:,1),xyseg(:,2),xy2,obs);
d = cumsum([0;sqrt(sum(diff(xyseg).^2,2))]);
dis = trapz(d,zseg);

%% finding risk factor
% close all
% for j = 1:length(obs.x_obs)
% rectangle('Position',[obs.x_obs(j) obs.y_obs(j) obs.w_obs(j) obs.h_obs(j)]);
% end
% hold
% 
% scatter(xyseg(:,1),xyseg(:,2));
% Assuming uniform
risk = 1;
for i = 1:length(obs.x_obs)
    p(i,:) = unifcdf(obs.x_obs(i) - xyseg(:,1),-obs.rx/2,obs.rx/2,'upper').* ...
            unifcdf(xyseg(:,1) - obs.x_obs(i) - obs.w_obs(i),-obs.rx/2,obs.rx/2,'upper').*...
            unifcdf(obs.y_obs(i) - xyseg(:,2),-obs.ry/2,obs.ry/2,'upper').*...
            unifcdf( xyseg(:,2) - obs.y_obs(i) - obs.h_obs(i),-obs.ry/2,obs.ry/2,'upper');
end
risk = p;

return
