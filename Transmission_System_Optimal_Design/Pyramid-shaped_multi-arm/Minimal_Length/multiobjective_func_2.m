function F = multiobjective_func_2(x)
%x(1);      %ODs_1
%x(2);      %ODs_2
%x(3);      %lb
%x(4);      %alpha
%x(5);      %le

%x(6);      %beta_1
%x(7);      %beta_2

% Objectiv 1 - maximum twist angle of SS tube
F1 = -(1007958012753983*((5483660254958579*x(3))/576460752303423488 + (5483660254958579*x(5))/576460752303423488 + 38385621784710053/28823037615171174400))/(39582418599936000000000*pi*((x(2) + 1/5000)^4 - x(1)^4));

% Objectiv 2 - SS tube length
F2 = x(3) + x(5);         

F = [F1,F2];
end