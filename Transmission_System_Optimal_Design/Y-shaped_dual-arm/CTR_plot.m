%% 绘制优化后的管子设计结果
R = 1000*sol.lb/sol.alpha    %圆弧段半径， mm
P1 = [1000*lc,1000*a/2]
P2 = [1000*lc+R*sin(sol.alpha),1000*a/2+R*(1-cos(sol.alpha))]
P3 = [1000*lc+R*sin(sol.alpha)+1000*(le_min+delta_L_1)*cos(sol.alpha),1000*a/2+R*(1-cos(sol.alpha))+1000*(le_min+delta_L_2)*sin(sol.alpha)]
P4 = [1000*lc+R*sin(sol.alpha)+1000*(le_min+delta_L_2)*cos(sol.alpha),1000*a/2+R*(1-cos(sol.alpha))+1000*(le_min+delta_L_2)*sin(sol.alpha)]

% x1 = 0:0.1:P1(1);
% y1 = P1(2);
% plot(x1,y1)
% hold on
% x2 = P1(1):0.1:P2(1);
% y2 = 3.6+R-sqrt(R^2-(x2-80)^2);
% plot(x2,y2)

x=0:0.1:P4(1);
for i=1:length(x)
    if i<length(0:0.1:P1(1))
        y(i) = P1(2);
    end
    
    if (i>=length(0:0.1:P1(1))) && (i<length(0:0.1:P2(1)))
        y(i) = P1(2)+R-sqrt(R^2-(x(i)-P1(1))^2);
    end 
    
    if (i>=length(0:0.1:P2(1))) && (i<=length(0:0.1:P4(1)))
        y(i) = tan(sol.alpha)*(x(i)-P2(1))+P2(2);
    end
end
plot(x,y);
hold on
plot(x,-y);
hold on
line([P1(1) P1(1)],[-35 35],'linestyle','--','Color','b');
line([P2(1) P2(1)],[-35 35],'linestyle','--','Color','b');
line([P3(1) P3(1)],[-35 35],'linestyle','--','Color','b');
axis equal;
xlim([0 500]);
ylim([-35 35]);
xlabel('x/mm');
ylabel('y/mm');

