function Xdot = model(t,in2,in3)
%MODEL
%    Xdot = MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 22:37:47

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th7 = in3(4,:);
th8 = in3(5,:);
th11 = in3(6,:);
th12 = in3(7,:);
th13 = in3(8,:);
th14 = in3(9,:);
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
t5 = x7+3.67890485041078e+1;
t6 = x8+4.9999999406895e+1;
t7 = x4+4.9999392660528e+1;
t8 = x6+1.68751497203814e+1;
t9 = x9+3.5997576222139;
t10 = x5+5.84877021264063;
t14 = x3+3.53407066939944e-1;
t24 = x2+1.00593104953139e-4;
t11 = 1.0./t5;
t12 = 1.0./t6;
t13 = 1.0./t7;
t15 = 1.0./t8;
t18 = 1.0./t10;
t19 = 1.0./t9;
t22 = 1.0./t14;
t27 = 1.0./t24;
t16 = t12.*x8;
t17 = t13.*x4;
t20 = t11.*x7;
t21 = t15.*x6;
t23 = t18.*x5;
t25 = t19.*x9;
t26 = t22.*x3;
t28 = t27.*x2;
t29 = t16+t17+t20+t21+t23+t25+t26+t28;
mt1 = [-x1.*(th2-(t29.*th1)./8.0);-x1.*(t2+(t29.*th1)./(th3.*8.0));-x1.*(t3+t29.*th1.*1.049082969790315e-1);-x1.*(t4+t29.*th1.*2.500000059310496e-3);x1.*(t2+t3+t4-t29.*th1.*2.680384886305695e-1-th14.*x5);x1.*(th14.*x6-(t29.*th1)./(th7.*8.0));x1.*(th14.*x7-(t29.*th1)./(th8.*8.0))];
mt2 = [x1.*(t4-t29.*th1.*2.500000041938851e-3);x1.*(t3-t29.*th1.*2.105773257114671e-1);t2.*x1];
Xdot = [mt1;mt2];