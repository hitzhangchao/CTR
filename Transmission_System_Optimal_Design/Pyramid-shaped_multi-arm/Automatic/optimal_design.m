%% Optimal Design Framework for Pyramid-shaped Transmission System
% 主优化程序
% by Zhang Chao
% Date：2022/10/20

clear;clc;
tic
%% Load CTR params
CTR_params;

%% Formulation of the multiobjective optimization problem
func = @objective_func;                 %objective function
nvars = 5;                              %numbers of optimization variables
Aineq = [-1,1,0,0,0];                   %linear constraint: Ax<=b
bineq = -2*w-2*delta;                   %linear inequally constraint: Ax<=b
Aeq = [];                               %linear equally constraints:Aeq x = beq
beq = [];
lbd = [ODn_1+2*b+2*w,max(ODn_1+2*b-2*delta,ODn_2+2*b+2*delta),0.1,deg2rad(3),le_min];    % lower bound
ubd = [min(a,D-a/sin(pi/n)),min(a,D-a/sin(pi/n))-2*w-2*delta,0.5,deg2rad(10),5*le_min];                 % upper bound
nonlcon = @design_constraints;               %nonlinear constraints
PopulationSize_Data= 1000;                      %gamultiobj setting parameters
FunctionTolerance_Data = 1e-9;                  %gamultiobj setting parameters
ConstraintTolerance_Data = 1e-9;                %gamultiobj setting parameters

options = optimoptions('gamultiobj');           %default options
%modify options setting
options = optimoptions(options,'PopulationSize', PopulationSize_Data);
options = optimoptions(options,'FunctionTolerance', FunctionTolerance_Data);
options = optimoptions(options,'ConstraintTolerance', ConstraintTolerance_Data);
options = optimoptions(options,'CreationFcn', @gacreationnonlinearfeasible);
options = optimoptions(options,'CrossoverFcn', {  @crossoverintermediate [] });
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @gaplotpareto });   % plot pareto front
toc
% run galmultiobj function and get the pareto front
[x_gal,fval_gal] = gamultiobj(func,nvars,Aineq,bineq,Aeq,beq,lbd,ubd,nonlcon,options);
