function [p,T_tip,d_tip,ss,x] = ctr_fk_compliant(len_ex,theta,len,len_cu,uy_star,I)
% ETR柔性臂正运动学计算函数
% 输入：
    % len_ex -- 管子伸出{0}部分数组 - 运动学输入
    % theta --  管子旋转角度数组 - 运动学输入
    % len --    管子总长数组
    % len_cu -- 管子预弯曲段长度数组
    % uy_star --预曲率y分量
    % I --      各管惯性矩
    % x --      各管theta_dot(0)数组
% 输出：
    % p --      CTR柔性臂形状位置点
    % T_tip --  CTR柔性臂末端位姿
    % d_tip --  segment后各管顶点弧长位置
    % ss --     各弧长计算点(返回出去用于绘图)


%% Solving BVP
x0 = zeros(length(theta),1);
%x0=[0;0;0]; 
A=[]; b=[];
% 优化变量x：为各管输入端theta_dot(s=0)
% 优化目标COST：各管末端theta_dot值
% 初始值x0：theta_dot(0)，优化过程会自动改变x0
% 思想：打靶法的思想，将边界值问题转化为初值问题
x = fmincon(@(x)COST(x,len_ex,theta,len,len_cu,uy_star,I),x0,A,b);      % x为各管theta_dot(0)数组


%% 计算形状和末端位姿
[p,R_end,ss,d_start,d_tip] = ctr_shape(x,len_ex,theta,len,len_cu,uy_star,I);

T_tip = [R_end,p(:,end); [0,0,0,1]];     % tip pose
%T_tip = zeros(4,4);
end





%% 优化目标 函数
function cost = COST(theta_dot_0,len_ex,theta,len,len_cu,uy_star,I)
% 优化变量：各管初始theta_dot(s=0)
% 优化目标：各管tip处theta_dot(s=Li)的平方（最小化末端theta_dot，打靶法的思想，将边界值问题转化为初值问题）
% INPUT:
    % theta_dot_0 -- 各管输入端theta_dot(s=0)
% OUTPUT:
    % cost -- 各管末端theta_dot(s=L)的平方

%CTR_params;                 % 加载物理参数
n = length(len_ex);          % 管子总数

%*** Segmenting ***%
[S,d_start,d_tip,uy_star,I] = segment(len,len_ex,len_cu,uy_star,I);
m = length(S);              % 分段段数

SS = zeros(1,m);
for i=1:m
    SS(i) = sum(S(1:i));    % S是各段的长度，SS是从最内管proximal点起到各段的累计长度
end

% 只取{0}之后的segments
SS_ex = SS(SS+min(d_start)>0) + min(d_start);   % SS_ex表示{0}开始，到各段末端的弧长位置
m_ex = length(SS_ex);       % 伸出{0}的段数
temp = zeros(n,m_ex);       % arm的管数量在ETR中最好分别定义一下
uy_star_ex = temp; I_ex = temp; % 预定义数组
for i=1:n
    uy_star_ex(i,:) = uy_star(i,SS+min(d_start)>0);
    I_ex(i,:) = I(i,SS+min(d_start)>0);
end

%*** Solving ODE ***%
span = [0 SS_ex];                       % 从{0}到各段末端弧长位置（包括0）
ss=[];
theta_dot=[];                       
theta = theta-d_start.*theta_dot_0';    % {0}处各管的转动角度（这样就统一了输入端到0位置，考虑了0位置之前管的扭转）

% 循环求解各segment的ODE
for i=1:m_ex
    s_span = [span(i),span(i+1)];       % 某段的弧长位置区间
    y0 = zeros(2*n+12,1);               % 微分变量y初始值
    y0(1:n) = theta;                    % 各管转动角度theta
    y0(n+1:2*n) = theta_dot_0;          % 各管初始theta_dot
    y0(2*n+1:2*n+12) = zeros(12,1);     % 优化目标中用不到p和R，故可简化（直接令其为0）
    
    %每个segment区间以y0为初始值计算一次CTR微分方程ODE,返回计算弧长位置s和微分变量y
    
    [s,y] = ode23(@(s,y) ctr_ode(s,y,uy_star_ex(:,i),I_ex(:,i)),s_span,y0);
    ss = [ss;s];                        % 不断加上新求出的各segment的弧长计算点
    theta_dot = [theta_dot;y(:,n+1:2*n)];
    theta = y(end,1:n)';                % 当前segment末端theta值，作为下一segment初值
    theta_dot_0 = y(end,n+1:2*n)';      % 将当前segment末端theta_dot值作为下一segment初值 
end

% 找出各管tip处theta_dot值
theta_dot_tip = zeros(n,1);
for i=1:n
    [~,index] = min(abs(ss-d_tip(i)));  % 从各计算弧长点ss找出各管顶点是哪一个index
    theta_dot_tip(i) = theta_dot(index,i);
end

cost = theta_dot_tip'*theta_dot_tip;    % 优化目标函数 - 各管tip处theta_dot值的平方
end



%% CTR的ODE 函数
function dydt = ctr_ode(~,y,uy_star_ex,I_ex)
% CTR微分方程（自己推导的形式）
% INPUT:
    % y --  各变量数组
        % y(1:n) -- 各管转动角度theta
        % y(n+1:2*n) -- 各管theta_dot
        % y(2*n+1:2*n+3) -- 最内管位置
        % y(2*n+4:2*n+12) -- 最内管旋转矩阵各元素
    % uy_star_ex -- 伸出段预曲率y分量
    % I_ex --       伸出段惯性矩
% OUTPUT:
    % dydt -- 一阶微分方程组
    
%=== 减少调用时间 - CTR_parms用到的变量直接写出 ===%
% CTR_params;               % 加载物理参数
n = size(uy_star_ex,1);     % number of tubes
E = 60E9;                   % NiTi弹性模量，Pa
poisson_rate = 0.3;         % NiTi泊松比

dydt = zeros(length(y),1);
Kb_ex = I_ex*E;     % 伸出各段弯曲刚度(如果要标定各管的E, 可在这里修改)

% 表示dydt(1:n)
dydt(1:n) = y(n+1:2*n);

% 表示dydt(n+1:2*n)
for i=1:n
    c_temp = (1+poisson_rate)*uy_star_ex(i)/sum(Kb_ex);
    sum_temp = 0;
    for j=1:n
        sum_temp = sum_temp + Kb_ex(j)*uy_star_ex(j)*sin(y(i)-y(j));
    end
    dydt(n+i) = c_temp*sum_temp;
end

% 表示dydt(2*n+1:2*n+3)
R = reshape(y(2*n+4:2*n+12),[3,3])';    % 最内管旋转矩阵，将1x9旋转矩阵元素还原为3x3矩阵
p_dot = R*[0;0;1];                      % p_dot = R*e3
dydt(2*n+1:2*n+3) = p_dot;

% 表示dydt(2*n+4:2*n+12)
ux = -sum(Kb_ex.*uy_star_ex.*sin(y(1:n)))/sum(Kb_ex);   % Centerline曲率分量
uy = sum(Kb_ex.*uy_star_ex.*cos(y(1:n)))/sum(Kb_ex);
un_x = ux*cos(y(n)) + uy*sin(y(n));     % 最内管曲率分量（非centerline）
un_y = -ux*sin(y(n)) + uy*cos(y(n));
un_z = y(2*n);

un_hat = skew([un_x;un_y;un_z]);       
R_dot = R*un_hat;                % R_dot = R*u_hat
R_dot = reshape(R_dot',[9,1]);  % 将3x3的R_dot重构为9x1的向量
dydt(2*n+4:2*n+12) = R_dot;
end



%% CTR形状和末端位姿 函数
function [p,R_end,ss,d_start,d_tip] = ctr_shape(theta_dot_0,len_ex,theta,len,len_cu,uy_star,I)
% 带入求解ODE后得到的theta_dot(s=0),求解CTR形状和末端位姿
% INPUT:
    % q --  运动学输入
    % theta_dot_0 --    输入端theta_dot(s=0)
% OUTPUT:
    % p --  最内管空间位置
    % R_end --  最内管末端旋转矩阵
    % ss -- 计算的弧长位置

n = length(len);            % 管子总数
theta_n = theta(end);       % 最内管初始旋转角度

%*** Segmenting ***%
[S,d_start,d_tip,uy_star,I] = segment(len,len_ex,len_cu,uy_star,I);
m = length(S);              % 分段段数

SS = zeros(1,m);
for i=1:m
    SS(i) = sum(S(1:i));    % S是各段的长度，SS是从最内管proximal点起到各段的累计长度
end

% 只取{0}之后的segments
SS_ex = SS(SS+min(d_start)>0) + min(d_start);   % SS_ex表示{0}开始，到各段末端的弧长位置
m_ex = length(SS_ex);        % 伸出{0}的段数
temp = zeros(n,m_ex); 
uy_star_ex = temp; I_ex = temp; % 预定义数组
for i=1:n
    uy_star_ex(i,:) = uy_star(i,SS+min(d_start)>0);
    I_ex(i,:) = I(i,SS+min(d_start)>0);
end

%*** Solving ODE ***%
span = [0 SS_ex];                   % 从{0}到各段末端弧长位置（包括0）
ss=[];p=[];
%theta_s=[];
%theta_dot=[];
p0 = [0;0;0];
R0 = rotz(theta_n);
R0 = reshape(R0',[9,1]);                % 将3x3的R0重构为9x1的向量
theta = theta-d_start.*theta_dot_0';    % {0}处各管的转动角度（这样就统一了输入端到0位置）

% 循环求解各segment的ODE
for i=1:m_ex
    s_span = [span(i),span(i+1)];       % 某段的弧长位置区间
    y0 = zeros(2*n+12,1);               % 微分变量y初始值
    y0(1:n) = theta;                    % 各管转动角度theta
    y0(n+1:2*n) = theta_dot_0;          % 各管初始theta_dot
    y0(2*n+1:2*n+12) = [p0;R0];         % 各管位置和旋转矩阵
    
    %每个segment区间以y0为初始值计算一次CTR微分方程ODE,返回计算弧长位置s和微分变量y
    [s,y] = ode23(@(s,y) ctr_ode(s,y,uy_star_ex(:,i),I_ex(:,i)),s_span,y0);
    
    ss = [ss;s];                            % 不断加上新求出的各segment的弧长计算点
    %theta_s = [theta_s;y(:,1:n)];           % 不断加上新求出的各segment弧长点对应theta
    %theta_dot = [theta_dot;y(:,n+1:2*n)];   % 不断加上新求出的各segment弧长点对应theta_dot
    p = [p;y(:,2*n+1:2*n+3)];               % 不断加上新求出的各segment最内管位空间位置点p
    
    theta = y(end,1:n)';                    % 当前segment末端theta值，作为下一segment初值
    theta_dot_0 = y(end,n+1:2*n)';          % 当前segment末端theta_dot值作为下一segment初值 
    p0 = y(end,2*n+1:2*n+3)';               % 末端p0
    R0 = y(end,2*n+4:2*n+12)';              % 末端R0
end

R_end = reshape(R0,[3,3])';
p = p'; 
end
