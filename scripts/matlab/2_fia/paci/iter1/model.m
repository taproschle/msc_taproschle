function Xdot = model(t,in2,in3)
%MODEL
%    Xdot = MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 23:13:55

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th4 = in3(4,:);
th5 = in3(5,:);
th6 = in3(6,:);
th7 = in3(7,:);
th8 = in3(8,:);
th9 = in3(9,:);
th10 = in3(10,:);
th11 = in3(11,:);
th12 = in3(12,:);
th13 = in3(13,:);
th14 = in3(14,:);
th18 = in3(15,:);
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
t5 = th18+x5;
t6 = x4+5.0e+1;
t7 = x8+5.0e+1;
t14 = x3+4.60072948671266e+1;
t15 = x6+3.92916983687822e+1;
t16 = x9+2.58902112996201e+1;
t17 = x2+3.27288990615856e+1;
t18 = x7+4.87007202410084e+1;
t8 = 1.0./t5;
t9 = 1.0./t6;
t10 = 1.0./t7;
t19 = 1.0./t14;
t20 = 1.0./t15;
t21 = 1.0./t16;
t22 = 1.0./t17;
t23 = 1.0./t18;
t11 = t8.*x5;
t12 = t9.*x4;
t13 = t10.*x8;
t24 = t19.*x3;
t25 = t20.*x6;
t26 = t21.*x9;
t27 = t22.*x2;
t28 = t23.*x7;
t29 = t11+t12+t13+t24+t25+t26+t27+t28;
Xdot = [-x1.*(th2-(t29.*th1)./8.0);-x1.*(t2+(t29.*th1)./(th3.*8.0));-x1.*(t3+(t29.*th1)./(th4.*8.0));-x1.*(t4+(t29.*th1)./(th5.*8.0));x1.*(t2+t3+t4-th14.*x5-(t29.*th1)./(th6.*8.0));x1.*(th14.*x6-(t29.*th1)./(th7.*8.0));x1.*(th14.*x7-(t29.*th1)./(th8.*8.0));x1.*(t4-(t29.*th1)./(th9.*8.0));x1.*(t3-(t29.*th1)./(th10.*8.0));t2.*x1];