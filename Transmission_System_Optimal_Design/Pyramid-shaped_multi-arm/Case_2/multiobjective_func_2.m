function F = multiobjective_func_2(x)
%x(1);      %ODs_1
%x(2);      %ODs_2
%x(3);      %ODs_3
%x(4);      %lb
%x(5);      %alpha
%x(6);      %le

%x(7);      %beta_1
%x(8);      %beta_2
%x(9);      %beta_3

% Objectiv 1 - maximum twist angle of SS tube
F1 = (1007958012753983*((93*x(4))/10000 + (93*x(6))/10000 + 93/40000))/(39582418599936000000000*pi*(x(3)^4 - 2812425566334057/618970019642690137449562112));

% Objectiv 2 - SS tube length
F2 = x(4) + x(6) + 0.220;         

F = [F1,F2];
end