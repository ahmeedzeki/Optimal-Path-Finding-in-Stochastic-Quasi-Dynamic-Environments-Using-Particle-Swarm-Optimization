clear
close all
clc
rng(3)
%% Defining global parameters
x_des = 45; % the target (destination) x axis point
y_des = 45; % the target (destination) y axis point
x0 = 5;
y0 = 5;

MaxIt = 250;
w=2;            % Inertia Weight
wdamp=.99;     % Inertia Weight Damping Ratio
c1=.3;         % Personal Learning Coefficient
c2=.2;         % Global Learning Coefficient
c3 = 0;         % Personal Learning Coefficient (worst)
VelMin = -5;
VelMax = 5;
VarMax = 50;
VarMin = 0;

render = 0; % render=1 update the figure
%% Construct the obs
obs_diff = 1; %1:simple scenario, 2:moderately difficult, 3: difficult scenario
switch obs_diff
    case 1
       
        obs.x_obs = [ 28 15];
        obs.y_obs = [ 25 5];
        obs.w_obs = [ 25 9];
        obs.h_obs = [ 13 25];
        obs.rx = 3;
        obs.ry = 3;
        scenario = 'simple';
      
    case 2
        obs.x_obs = [15  15];
        obs.y_obs = [35  5];
        obs.w_obs = [5  9];
        obs.h_obs = [20  25];
        obs.rx = 3;
        obs.ry = 3;
        scenario = 'moderate';
    case 3
        obs.x_obs = [15 28 15];
        obs.y_obs = [35 25 5];
        obs.w_obs = [5 25 9];
        obs.h_obs = [20 13 25];
        obs.rx = 3;
        obs.ry = 3;
        scenario = 'hard';
       
    otherwise
        disp('Error ')
        return
end

    
    



%% Defining pso parameters
nPop = 500;
alpha = 10000;
NP = 2; % Number of points

VarSize = 2;

% Initialization
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.risk_factor = [];
empty_particle.Worst.Cost = [];
particle=repmat(empty_particle,nPop,1);
GlobalBest.Cost=inf;
for i=1:nPop
    
    % Initialize Position
    particle(i).Position = [50*rand(NP,1) 50*rand(NP,1)];
    % Initalize worst positon
    
    % Initialize Velocity
    particle(i).Velocity=.5-rand(1,VarSize);
    
    % Evaluation
    Vp = [[x0;y0] particle(i).Position' [x_des;y_des]];
    cst = 0;
    rft = 0;
    for lm = 1:(length(Vp)-1)
        [~,cs_temp,rf] = line_integral2_with_risk_factor(Vp(:,lm)',Vp(:,lm+1)',obs);
        cst = cst + cs_temp;
        rft = rft + sum(rf(:));
    end
    particle(i).Cost = cst+alpha*rft;
    particle(i).risk_factor = rft;
    % worst cost
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end

BestCost=zeros(MaxIt,1);

%% Creating the scene
x_start = 0; % starting position of x axis scene area
x_end = 50; % ending position of x axis scene area
y_start = 0; % starting position of y axis scene area
y_end = 50; % ending position of y axis scene area
rec_width = 2;

figure
GlobalBest.index = 1;
for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        velocity = w*particle(i).Velocity ...
            +c1*rand(1,VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(1,VarSize).*(GlobalBest.Position-particle(i).Position);
        particle(i).Velocity =  velocity ;
%         display(particle(i).Velocity)
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        new_pos = particle(i).Position + particle(i).Velocity;
            particle(i).Position = new_pos;
            Vp = [[x0;y0] particle(i).Position' [x_des;y_des]];
            cst = 0;
            rft = 0;
            for lm = 1:(length(Vp)-1)
                [obs_t,cs_temp,rf] = line_integral2_with_risk_factor(Vp(:,lm)',Vp(:,lm+1)',obs);
                cst = cst + cs_temp;
                rft = rft + sum(rf(:));
                
            end
            particle(i).Cost = cst+alpha*rft;
            particle(i).risk_factor = rft;
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        
        
        
        % Update Personal Best
         if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                GlobalBest.index = i;
                
            end
            
        end
        
 
       
        
    end
    
    BestCost(it)=GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(particle(GlobalBest.index).Best.Cost) ': Risk Factor = ' num2str(particle(GlobalBest.index).risk_factor)]);
    
    w=w*wdamp;
       
    % Updating figure
    if render
        clf()
        hold
        xax = [x0 particle(GlobalBest.index).Best.Position(:,1)' x_des];
        yax = [y0 particle(GlobalBest.index).Best.Position(:,2)' y_des];
        plot(xax,yax);

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
        pause(.1)
        axis equal

        hold
    end
end
BestSol = GlobalBest;

    

%%
clf()
    hold
        xax = [x0 particle(GlobalBest.index).Best.Position(:,1)' x_des];
        yax = [y0 particle(GlobalBest.index).Best.Position(:,2)' y_des];
       plot(xax,yax);
   
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
    file_path_fig = ['results_alpha_',num2str(alpha),'_np_',num2str(NP),'\\','npop_',num2str(nPop)];
    file_path_ws = ['results_alpha_',num2str(alpha),'_np_',num2str(NP),'\\','npop_',num2str(nPop)];
    file_path_txt = ['results_alpha_',num2str(alpha),'_np_',num2str(NP),'\\','results.txt'];
    file_convergence = ['results_alpha_',num2str(alpha),'_np_',num2str(NP),'\\','results_convergence_',scenario,'.txt'];
    savefig(file_path_fig)
    save(file_path_ws)
    
fid = fopen(file_path_txt, 'a+');
v = ['n_pop = ',num2str(nPop),'\t',' Best Cost = ' num2str(particle(GlobalBest.index).Best.Cost) ': Risk Factor = ' num2str(particle(GlobalBest.index).risk_factor),'\n'];
fprintf(fid,v);
fclose(fid);

fid = fopen(file_convergence, 'a+');
for i = 1:MaxIt
    v = [ 'Np = ',num2str(NP),'\t','alpha = ',num2str(alpha),'\t','iteration = ',num2str(i),'\t',' Best Cost = ' num2str(BestCost(i)),'\n'];
    fprintf(fid,v);
end
fclose(fid);