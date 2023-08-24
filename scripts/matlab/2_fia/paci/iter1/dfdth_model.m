function dfdthSym = dfdth_model(t,in2,in3)
%DFDTH_MODEL
%    dfdthSym = DFDTH_MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 23:13:56

th1 = in3(1,:);
th3 = in3(3,:);
th4 = in3(4,:);
th5 = in3(5,:);
th6 = in3(6,:);
th7 = in3(7,:);
th8 = in3(8,:);
th9 = in3(9,:);
th10 = in3(10,:);
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
t2 = th18+x5;
t3 = x1.*x2;
t4 = x1.*x3;
t5 = x1.*x4;
t6 = 1.0./th3;
t7 = 1.0./th4;
t8 = 1.0./th5;
t9 = 1.0./th6;
t10 = 1.0./th7;
t11 = 1.0./th8;
t12 = 1.0./th9;
t13 = 1.0./th10;
t14 = x4+5.0e+1;
t15 = x8+5.0e+1;
t23 = x3+4.60072948671266e+1;
t24 = x6+3.92916983687822e+1;
t25 = x9+2.58902112996201e+1;
t26 = x2+3.27288990615856e+1;
t27 = x7+4.87007202410084e+1;
t16 = 1.0./t2;
t18 = 1.0./t14;
t19 = 1.0./t15;
t28 = 1.0./t23;
t29 = 1.0./t24;
t30 = 1.0./t25;
t31 = 1.0./t26;
t32 = 1.0./t27;
t17 = t16.^2;
t20 = t16.*x5;
t21 = t18.*x4;
t22 = t19.*x8;
t33 = t28.*x3;
t34 = t29.*x6;
t35 = t30.*x9;
t36 = t31.*x2;
t37 = t32.*x7;
t38 = t20+t21+t22+t33+t34+t35+t36+t37;
mt1 = [x1.*(t20./8.0+t21./8.0+t22./8.0+t33./8.0+t34./8.0+t35./8.0+t36./8.0+t37./8.0),t6.*t38.*x1.*(-1.0./8.0),t7.*t38.*x1.*(-1.0./8.0),t8.*t38.*x1.*(-1.0./8.0),t9.*t38.*x1.*(-1.0./8.0),t10.*t38.*x1.*(-1.0./8.0),t11.*t38.*x1.*(-1.0./8.0),t12.*t38.*x1.*(-1.0./8.0),t13.*t38.*x1.*(-1.0./8.0),0.0,-x1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t6.^2.*t38.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t7.^2.*t38.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t8.^2.*t38.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
mt2 = [0.0,(t9.^2.*t38.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t10.^2.*t38.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t11.^2.*t38.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t12.^2.*t38.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t13.^2.*t38.*th1.*x1)./8.0,0.0,0.0,-t3,0.0,0.0,t3,0.0,0.0,0.0,0.0,t3,0.0,0.0,-t4,0.0,t4,0.0,0.0,0.0,t4,0.0,0.0,0.0,0.0,-t5,t5,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,0.0,-x1.*x5,x1.*x6,x1.*x7,0.0,0.0,0.0,t17.*th1.*x1.*x5.*(-1.0./8.0),(t6.*t17.*th1.*x1.*x5)./8.0,(t7.*t17.*th1.*x1.*x5)./8.0,(t8.*t17.*th1.*x1.*x5)./8.0];
mt3 = [(t9.*t17.*th1.*x1.*x5)./8.0,(t10.*t17.*th1.*x1.*x5)./8.0,(t11.*t17.*th1.*x1.*x5)./8.0,(t12.*t17.*th1.*x1.*x5)./8.0,(t13.*t17.*th1.*x1.*x5)./8.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3],10,15);