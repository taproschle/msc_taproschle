function dfdthSym = dfdth_model(t,in2,in3)
%DFDTH_MODEL
%    dfdthSym = DFDTH_MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 23:16:08

th1 = in3(1,:);
th3 = in3(3,:);
th4 = in3(4,:);
th6 = in3(5,:);
th7 = in3(6,:);
th8 = in3(7,:);
th18 = in3(12,:);
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
t8 = 1.0./th6;
t9 = 1.0./th7;
t10 = 1.0./th8;
t11 = x4+5.0e+1;
t12 = x8+5.0e+1;
t24 = x3+4.60072948671266e+1;
t25 = x6+3.92916983687822e+1;
t26 = x9+2.58902112996201e+1;
t27 = x2+3.27288990615856e+1;
t28 = x7+4.87007202410084e+1;
t13 = 1.0./t2;
t15 = 1.0./t11;
t16 = 1.0./t12;
t29 = 1.0./t24;
t30 = 1.0./t25;
t31 = 1.0./t26;
t32 = 1.0./t27;
t33 = 1.0./t28;
t14 = t13.^2;
t17 = t13.*x5;
t18 = t15.*x4;
t19 = t16.*x8;
t34 = t29.*x3;
t35 = t30.*x6;
t36 = t31.*x9;
t37 = t32.*x2;
t38 = t33.*x7;
t20 = t17./4.0e+2;
t21 = t18./4.0e+2;
t22 = t19./4.0e+2;
t23 = (t14.*th1.*x1.*x5)./4.0e+2;
t39 = t34./4.0e+2;
t40 = t35./4.0e+2;
t41 = t36./4.0e+2;
t42 = t37./4.0e+2;
t43 = t38./4.0e+2;
t44 = t17+t18+t19+t34+t35+t36+t37+t38;
t45 = t20+t21+t22+t39+t40+t41+t42+t43;
t46 = t45.*x1;
t47 = -t46;
et1 = t17.*2.500391799741333e-3+t18.*2.500391799741333e-3+t19.*2.500391799741333e-3+t34.*2.500391799741333e-3+t35.*2.500391799741333e-3+t36.*2.500391799741333e-3;
et2 = t37.*2.500391799741333e-3+t38.*2.500391799741333e-3;
mt1 = [x1.*(t17./8.0+t18./8.0+t19./8.0+t34./8.0+t35./8.0+t36./8.0+t37./8.0+t38./8.0),t6.*t44.*x1.*(-1.0./8.0),t7.*t44.*x1.*(-1.0./8.0),t47,t8.*t44.*x1.*(-1.0./8.0),t9.*t44.*x1.*(-1.0./8.0),t10.*t44.*x1.*(-1.0./8.0),t47,-x1.*(et1+et2),0.0,-x1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t6.^2.*t44.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t7.^2.*t44.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t8.^2.*t44.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t9.^2.*t44.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0];
mt2 = [0.0,0.0,0.0,0.0,0.0,(t10.^2.*t44.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,-t3,0.0,0.0,t3,0.0,0.0,0.0,0.0,t3,0.0,0.0,-t4,0.0,t4,0.0,0.0,0.0,t4,0.0,0.0,0.0,0.0,-t5,t5,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,0.0,-x1.*x5,x1.*x6,x1.*x7,0.0,0.0,0.0,t14.*th1.*x1.*x5.*(-1.0./8.0),(t6.*t14.*th1.*x1.*x5)./8.0,(t7.*t14.*th1.*x1.*x5)./8.0,t23,(t8.*t14.*th1.*x1.*x5)./8.0,(t9.*t14.*th1.*x1.*x5)./8.0,(t10.*t14.*th1.*x1.*x5)./8.0,t23,t14.*th1.*x1.*x5.*2.500391799741333e-3,0.0];
dfdthSym = reshape([mt1,mt2],10,12);
