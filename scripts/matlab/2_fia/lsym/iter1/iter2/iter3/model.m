function Xdot = model(t,in2,in3)
%MODEL
%    Xdot = MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 23:09:30

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th5 = in3(4,:);
th11 = in3(5,:);
th12 = in3(6,:);
th13 = in3(7,:);
th14 = in3(8,:);
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
t5 = x2+5.0e+1;
t6 = x4+5.0e+1;
t11 = x6+4.14964190890544e+1;
t12 = x3+2.19491164389871e+1;
t13 = x7+4.55219841529752e+1;
t14 = x8+4.99998982430964e+1;
t15 = x9+4.89808949505372e+1;
t16 = x5+2.17266488267788;
t7 = 1.0./t5;
t8 = 1.0./t6;
t17 = 1.0./t11;
t18 = 1.0./t15;
t19 = 1.0./t12;
t20 = 1.0./t13;
t21 = 1.0./t14;
t22 = 1.0./t16;
t9 = t7.*x2;
t10 = t8.*x4;
t23 = t17.*x6;
t24 = t18.*x9;
t25 = t19.*x3;
t26 = t20.*x7;
t27 = t21.*x8;
t28 = t22.*x5;
t29 = t9+t10+t23+t24+t25+t26+t27+t28;
mt1 = [-x1.*(th2-(t29.*th1)./8.0);-x1.*(t2+(t29.*th1)./(th3.*8.0));-x1.*(t3+t29.*th1.*3.873074066115075e-1);-x1.*(t4+(t29.*th1)./(th5.*8.0));x1.*(t2+t3+t4-t29.*th1.*1.574567667912253e-1-th14.*x5);-x1.*(t29.*th1.*2.545016783061344e-3-th14.*x6)];
mt2 = [-x1.*(t29.*th1.*2.569220064269693e-3-th14.*x7);x1.*(t4-t29.*th1.*2.500000011131805e-3);x1.*(t3-t29.*th1.*4.91060599048407e-2);t2.*x1];
Xdot = [mt1;mt2];
