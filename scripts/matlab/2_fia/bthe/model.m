function Xdot = model(t,in2,in3)
%MODEL
%    Xdot = MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 22:40:44

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
th15 = in3(15,:);
th16 = in3(16,:);
th17 = in3(17,:);
th18 = in3(18,:);
th19 = in3(19,:);
th20 = in3(20,:);
th21 = in3(21,:);
th22 = in3(22,:);
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
t5 = th15+x2;
t6 = th16+x3;
t7 = th17+x4;
t8 = th18+x5;
t9 = th19+x6;
t10 = th20+x7;
t11 = th21+x8;
t12 = th22+x9;
t13 = 1.0./t5;
t14 = 1.0./t6;
t15 = 1.0./t7;
t16 = 1.0./t8;
t17 = 1.0./t9;
t18 = 1.0./t10;
t19 = 1.0./t11;
t20 = 1.0./t12;
t21 = t13.*x2;
t22 = t14.*x3;
t23 = t15.*x4;
t24 = t16.*x5;
t25 = t17.*x6;
t26 = t18.*x7;
t27 = t19.*x8;
t28 = t20.*x9;
t29 = t21+t22+t23+t24+t25+t26+t27+t28;
Xdot = [-x1.*(th2-(t29.*th1)./8.0);-x1.*(t2+(t29.*th1)./(th3.*8.0));-x1.*(t3+(t29.*th1)./(th4.*8.0));-x1.*(t4+(t29.*th1)./(th5.*8.0));x1.*(t2+t3+t4-th14.*x5-(t29.*th1)./(th6.*8.0));x1.*(th14.*x6-(t29.*th1)./(th7.*8.0));x1.*(th14.*x7-(t29.*th1)./(th8.*8.0));x1.*(t4-(t29.*th1)./(th9.*8.0));x1.*(t3-(t29.*th1)./(th10.*8.0));t2.*x1];
