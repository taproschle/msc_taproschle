function Xdot = model(t,in2,in3)
%MODEL
%    Xdot = MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 21:23:10

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th4 = in3(4,:);
th6 = in3(5,:);
th7 = in3(6,:);
th8 = in3(7,:);
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
t2 = th12.*x3;
t3 = th13.*x4;
t4 = x2+5.0e+1;
t5 = x4+5.0e+1;
t6 = x8+5.0e+1;
t13 = x2./1.0e+4;
t14 = x9+1.16185559907235;
t15 = x6+1.20089248843775e+1;
t16 = x7+1.86021006236492e+1;
t17 = x5+1.41124995168232;
t23 = x3+2.63462005905101e-3;
t7 = 1.0./t4;
t8 = 1.0./t5;
t9 = 1.0./t6;
t18 = 1.0./t14;
t19 = 1.0./t15;
t20 = 1.0./t16;
t21 = 1.0./t17;
t27 = 1.0./t23;
t10 = t7.*x2;
t11 = t8.*x4;
t12 = t9.*x8;
t22 = t18.*x9;
t24 = t19.*x6;
t25 = t20.*x7;
t26 = t21.*x5;
t28 = t27.*x3;
t29 = t10+t11+t12+t22+t24+t25+t26+t28;
t30 = (t29.*th1)./4.0e+2;
Xdot = [-x1.*(th2-(t29.*th1)./8.0);-x1.*(t13+(t29.*th1)./(th3.*8.0));-x1.*(t2+(t29.*th1)./(th4.*8.0));-x1.*(t3+t30);x1.*(t2+t3+t13-th14.*x5-(t29.*th1)./(th6.*8.0));x1.*(th14.*x6-(t29.*th1)./(th7.*8.0));x1.*(th14.*x7-(t29.*th1)./(th8.*8.0));x1.*(t3-t29.*th1.*2.501230076257898e-3);x1.*(t2-t30);t13.*x1];