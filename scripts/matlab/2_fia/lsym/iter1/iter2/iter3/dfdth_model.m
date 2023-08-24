function dfdthSym = dfdth_model(t,in2,in3)
%DFDTH_MODEL
%    dfdthSym = DFDTH_MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 23:09:31

th1 = in3(1,:);
th3 = in3(3,:);
th5 = in3(4,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
x6 = in2(6,:);
x7 = in2(7,:);
x8 = in2(8,:);
x9 = in2(9,:);
t2 = x1.*x2;
t3 = x1.*x3;
t4 = x1.*x4;
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
et1 = t9.*3.873074066115075e-1+t10.*3.873074066115075e-1+t23.*3.873074066115075e-1+t24.*3.873074066115075e-1+t25.*3.873074066115075e-1;
et2 = t26.*3.873074066115075e-1+t27.*3.873074066115075e-1+t28.*3.873074066115075e-1;
et3 = t9.*1.574567667912253e-1+t10.*1.574567667912253e-1+t23.*1.574567667912253e-1+t24.*1.574567667912253e-1+t25.*1.574567667912253e-1;
et4 = t26.*1.574567667912253e-1+t27.*1.574567667912253e-1+t28.*1.574567667912253e-1;
et5 = t9.*2.545016783061344e-3+t10.*2.545016783061344e-3+t23.*2.545016783061344e-3+t24.*2.545016783061344e-3+t25.*2.545016783061344e-3+t26.*2.545016783061344e-3;
et6 = t27.*2.545016783061344e-3+t28.*2.545016783061344e-3;
et7 = t9.*2.569220064269693e-3+t10.*2.569220064269693e-3+t23.*2.569220064269693e-3+t24.*2.569220064269693e-3+t25.*2.569220064269693e-3+t26.*2.569220064269693e-3;
et8 = t27.*2.569220064269693e-3+t28.*2.569220064269693e-3;
et9 = t9.*2.500000011131805e-3+t10.*2.500000011131805e-3+t23.*2.500000011131805e-3+t24.*2.500000011131805e-3+t25.*2.500000011131805e-3+t26.*2.500000011131805e-3;
et10 = t27.*2.500000011131805e-3+t28.*2.500000011131805e-3;
et11 = t9.*4.91060599048407e-2+t10.*4.91060599048407e-2+t23.*4.91060599048407e-2+t24.*4.91060599048407e-2+t25.*4.91060599048407e-2+t26.*4.91060599048407e-2;
et12 = t27.*4.91060599048407e-2+t28.*4.91060599048407e-2;
mt1 = [x1.*(t9./8.0+t10./8.0+t23./8.0+t24./8.0+t25./8.0+t26./8.0+t27./8.0+t28./8.0),(t29.*x1.*(-1.0./8.0))./th3,-x1.*(et1+et2),(t29.*x1.*(-1.0./8.0))./th5,-x1.*(et3+et4),-x1.*(et5+et6),-x1.*(et7+et8),-x1.*(et9+et10),-x1.*(et11+et12),0.0,-x1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t29.*th1.*1.0./th3.^2.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t29.*th1.*1.0./th5.^2.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t2,0.0,0.0,t2,0.0,0.0,0.0,0.0,t2,0.0,0.0,-t3,0.0,t3,0.0,0.0,0.0,t3,0.0,0.0];
mt2 = [0.0,0.0,-t4,t4,0.0,0.0,t4,0.0,0.0,0.0,0.0,0.0,0.0,-x1.*x5,x1.*x6,x1.*x7,0.0,0.0,0.0];
dfdthSym = reshape([mt1,mt2],10,8);
