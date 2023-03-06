function F = multiobjective_func(x)
%x(1);      %ODs_1
%x(2);      %ODs_2
%x(3);      %lb
%x(4);      %alpha
%x(5);      %le

%x(6);      %beta_1
%x(7);      %beta_2

% Objectiv 1 - maximum twist angle of SS tube
F1 = -(1007958012753983*((7215989436560773*x(3))/36028797018963968 + (7215989436560773*x(5))/36028797018963968 + 64943904929046957/1801439850948198400))/(39582418599936000000000*pi*((x(2) + 1/5000)^4 - x(1)^4));

% Objectiv 2 - SS tube length
F2 = x(3) + x(5);         

F = [F1,F2];
end