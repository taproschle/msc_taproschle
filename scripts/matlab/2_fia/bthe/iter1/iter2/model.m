function Xdot = model(t,in2,in3)
%MODEL
%    Xdot = MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 22:45:54

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th6 = in3(4,:);
th7 = in3(5,:);
th8 = in3(6,:);
th11 = in3(7,:);
th12 = in3(8,:);
th13 = in3(9,:);
th14 = in3(10,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
x6 = in2(6,:);
x7 = in2(7,:);
x8 = in2(8,:);
x9 = in2(9,:);
t2 = th11.*x2;
t3 = th12.*x3;
t4 = th13.*x4;
t5 = x8+3.7995870797099e+1;
t6 = x3+4.99999999906868e+1;
t7 = x6+2.19492689619356e+1;
t8 = x4+4.99999994988811e+1;
t9 = x9+1.86124538346878;
t10 = x7+1.66530861806818;
t11 = x5+3.84036420636605;
t14 = x2+8.16450007953423e-1;
t12 = 1.0./t6;
t13 = 1.0./t5;
t15 = 1.0./t7;
t16 = 1.0./t8;
t18 = 1.0./t11;
t20 = 1.0./t9;
t21 = 1.0./t10;
t24 = 1.0./t14;
t17 = t12.*x3;
t19 = t13.*x8;
t22 = t16.*x4;
t23 = t15.*x6;
t25 = t20.*x9;
t26 = t21.*x7;
t27 = t18.*x5;
t28 = t24.*x2;
t29 = t17+t19+t22+t23+t25+t26+t27+t28;
t30 = t29.*th1.*2.50000000046566e-3;
Xdot = [-x1.*(th2-(t29.*th1)./8.0);-x1.*(t2+(t29.*th1)./(th3.*8.0));-x1.*(t3+t29.*th1.*4.223924006036696e-2);-x1.*(t4+t30);x1.*(t2+t3+t4-th14.*x5-(t29.*th1)./(th6.*8.0));x1.*(th14.*x6-(t29.*th1)./(th7.*8.0));x1.*(th14.*x7-(t29.*th1)./(th8.*8.0));x1.*(t4-t30);x1.*(t3-t29.*th1.*2.593753788511239e-1);t2.*x1];