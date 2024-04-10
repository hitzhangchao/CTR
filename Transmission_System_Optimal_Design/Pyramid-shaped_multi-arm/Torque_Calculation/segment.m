function [S,d3,d1,uy_star,I] = segment(len,len_ex,len_cu,uy_star,I)
% 管子分段函数，从外到内标号逐渐增大（最外管编号为1）
% 输入：
    % len --    管子总长数组
    % len_ex -- 管子伸出{0}部分数组
    % len_cu -- 管子预弯曲段长度数组
    % uy_star --预曲率y分量
    % I --      各管惯性矩
% 输出：
    % S --      各段长度
    % d3 --     各管proximal段开始位置（为负值）
    % d1 --     各管tip弧长位置（即伸出长度）
    % uy_star --分段后的各段各管预曲率y分量
    % I --      分段后的各管惯性矩

k = length(len);          % 管子数量

d1 = len_ex;            % 各管tip弧长位置（即伸出长度）
d2 = d1 - len_cu;       % 各管直管段结束位置
d3 = len_ex - len;      % 各管proximal段开始位置（为负值）
points = [0 d3 d2 d1];  % 重要点的位置点

[L,index] = sort(points);       % 位置点排序，返回有序位置点L和点在原points数组中的index
S = 1E-6*floor(1E6*diff(L));    % 各segment长度（先相邻点位置用diff差分相减，增大数值用floor取个整，然后再缩小）
% used floor because diff command doesn't give absolute zero sometimes

% 增加proximal段管子位置判断，增强鲁棒性
for i=1:k-1
    if d3(i+1) >= d3(i)
        sprintf('inner tube is clashing into outer tubes')
        EE=zeros(k,length(S)); I=EE; uy_star=EE;
        return
    end
end

% 提前声明一些数组，对3管N segments情况，数值为3行N列
EE = zeros(k,length(S)); II=EE; UUy=EE;
for i=1:k
    a = find(index==i+1);       % find where tube begins（注意index是排序后的点在原points数组中的位置）
    b = find(index==k+i+1);     % find where tube curve starts
    c = find(index==2*k+i+1);   % find where tube ends

    if S(a)==0; a=a+1;  end
    if S(b)==0; b=b+1;  end
    if c<=length(S)
        if S(c)==0; c=c+1; end
    end

    II(i,a:c-1)=I(i);
    UUy(i,b:c-1)=uy_star(i);
end

len = S(~(S==0));  % get rid of zero lengthes
E=zeros(k,length(len)); I=E; uy_star=E;
 for i=1:k
    I(i,:)=II(i,~(S==0));
    uy_star(i,:)=UUy(i,~(S==0));
 end
 
 S = S(~(S==0));
end
