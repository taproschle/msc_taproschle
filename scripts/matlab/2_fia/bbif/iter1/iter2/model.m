function Xdot = model(t,in2,in3)
%MODEL
%    Xdot = MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 21:38:03

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th4 = in3(4,:);
th5 = in3(5,:);
th6 = in3(6,:);
th7 = in3(7,:);
th8 = in3(8,:);
th10 = in3(9,:);
th12 = in3(10,:);
th13 = in3(11,:);
th14 = in3(12,:);
th22 = in3(13,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
x6 = in2(6,:);
x7 = in2(7,:);
x8 = in2(8,:);
x9 = in2(9,:);
t2 = th12.*x3;
t3 = th13.*x4;
t4 = th22+x9;
t5 = x2+5.0e+1;
t10 = x2./1.0e+4;
t11 = x6+4.95506770873895e+1;
t12 = x5+4.37180387399094e+1;
t13 = x3+4.99954399431281e+1;
t14 = x8+1.88091400530181e+1;
t15 = x7+4.99999999939224e+1;
t21 = x4+1.23342243728317e-3;
t6 = 1.0./t4;
t7 = 1.0./t5;
t16 = 1.0./t14;
t17 = 1.0./t15;
t18 = 1.0./t11;
t19 = 1.0./t12;
t20 = 1.0./t13;
t27 = 1.0./t21;
t8 = t6.*x9;
t9 = t7.*x2;
t22 = t17.*x7;
t23 = t18.*x6;
t24 = t19.*x5;
t25 = t20.*x3;
t26 = t16.*x8;
t28 = t27.*x4;
t29 = t8+t9+t22+t23+t24+t25+t26+t28;
Xdot = [-x1.*(th2-(t29.*th1)./8.0);-x1.*(t10+(t29.*th1)./(th3.*8.0));-x1.*(t2+(t29.*th1)./(th4.*8.0));-x1.*(t3+(t29.*th1)./(th5.*8.0));x1.*(t2+t3+t10-th14.*x5-(t29.*th1)./(th6.*8.0));x1.*(th14.*x6-(t29.*th1)./(th7.*8.0));x1.*(th14.*x7-(t29.*th1)./(th8.*8.0));x1.*(t3-t29.*th1.*2.502801109457882e-3);x1.*(t2-(t29.*th1)./(th10.*8.0));t10.*x1];
