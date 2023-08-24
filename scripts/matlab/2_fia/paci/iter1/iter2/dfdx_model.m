function dfdxSym = dfdx_model(t,in2,in3)
%DFDX_MODEL
%    dfdxSym = DFDX_MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 23:16:07

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th4 = in3(4,:);
th6 = in3(5,:);
th7 = in3(6,:);
th8 = in3(7,:);
th11 = in3(8,:);
th12 = in3(9,:);
th13 = in3(10,:);
th14 = in3(11,:);
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
t2 = th11.*x2;
t3 = th12.*x3;
t4 = th13.*x4;
t5 = th18+x5;
t6 = 1.0./th3;
t7 = 1.0./th4;
t8 = 1.0./th6;
t9 = 1.0./th7;
t10 = 1.0./th8;
t11 = x4+5.0e+1;
t12 = x8+5.0e+1;
t28 = x3.*3.5184372088832e+13;
t29 = x2.*1.40737488355328e+14;
t30 = x6.*1.40737488355328e+14;
t31 = x7.*7.0368744177664e+13;
t32 = x9.*1.40737488355328e+14;
t35 = x3+4.60072948671266e+1;
t37 = x6+3.92916983687822e+1;
t38 = x9+2.58902112996201e+1;
t39 = x2+3.27288990615856e+1;
t40 = x7+4.87007202410084e+1;
t13 = 1.0./t5;
t15 = 1.0./t11;
t17 = 1.0./t12;
t36 = t28+1.618737781405593e+15;
t41 = t30+5.529814941637545e+15;
t42 = t32+3.643723311297265e+15;
t43 = t29+4.606183050562609e+15;
t44 = t31+3.427008523907503e+15;
t46 = 1.0./t35;
t51 = 1.0./t37;
t52 = 1.0./t38;
t53 = 1.0./t39;
t55 = 1.0./t40;
t14 = t13.^2;
t16 = t15.^2;
t18 = t17.^2;
t19 = t13.*x5;
t20 = t15.*x4;
t22 = t17.*x8;
t23 = -t15;
t45 = 1.0./t36.^2;
t47 = t46.^2;
t48 = 1.0./t42.^2;
t49 = 1.0./t43.^2;
t50 = 1.0./t44.^2;
t54 = t53.^2;
t56 = 1.0./t41.^2;
t57 = t46.*x3;
t59 = -t46;
t60 = t51.*x6;
t61 = t52.*x9;
t62 = t53.*x2;
t64 = t55.*x7;
t65 = -t53;
t21 = t16.*x4;
t24 = (t18.*th1.*x1)./8.0;
t25 = (t14.*th1.*th18.*x1)./4.0e+2;
t58 = t47.*x3;
t63 = t54.*x2;
t66 = t56.*th1.*x1.*1.945630664889582e+27;
t67 = t48.*th1.*x1.*1.28202116773434e+27;
t68 = t45.*th1.*x1.*1.42385681038062e+26;
t72 = t49.*th1.*x1.*1.620656583602661e+27;
t73 = t50.*th1.*x1.*6.028857152838025e+26;
t78 = t19+t20+t22+t57+t60+t61+t62+t64;
t26 = -t24;
t27 = -t25;
t33 = t21+t23;
t34 = th1.*(t15-t21).*(-1.0./4.0e+2);
t69 = -t66;
t70 = -t67;
t71 = -t68;
t74 = -t72;
t75 = -t73;
t76 = t58+t59;
t77 = t63+t65;
t79 = (t78.*th1)./4.0e+2;
t80 = -t79;
mt1 = [-th2+(t78.*th1)./8.0,-t2-(t6.*t78.*th1)./8.0,-t3-(t7.*t78.*th1)./8.0,-t4+t80,t2+t3+t4-th14.*x5-(t8.*t78.*th1)./8.0,th14.*x6-(t9.*t78.*th1)./8.0,th14.*x7-(t10.*t78.*th1)./8.0,t4+t80,t3-t78.*th1.*2.500391799741333e-3,t2,t49.*th1.*x1.*8.103282918013305e+28,-x1.*(th11+(t6.*th1.*(t53-t63))./8.0)];
mt2 = [t7.*t49.*th1.*x1.*-8.103282918013305e+28,t74,x1.*(th11-(t8.*th1.*(t53-t63))./8.0),t9.*t49.*th1.*x1.*-8.103282918013305e+28,t10.*t49.*th1.*x1.*-8.103282918013305e+28,t74,t49.*th1.*x1.*(-1.620910572734759e+27),th11.*x1];
mt3 = [t45.*th1.*x1.*7.119284051903098e+27,t6.*t45.*th1.*x1.*-7.119284051903098e+27,-x1.*(th12+(t7.*th1.*(t46-t58))./8.0),t71,x1.*(th12-(t8.*th1.*(t46-t58))./8.0),t9.*t45.*th1.*x1.*-7.119284051903098e+27,t10.*t45.*th1.*x1.*-7.119284051903098e+27,t71];
mt4 = [x1.*(th12-th1.*(t46-t58).*2.500391799741333e-3),0.0,t16.*th1.*x1.*(2.5e+1./4.0),t6.*t16.*th1.*x1.*(-2.5e+1./4.0),t7.*t16.*th1.*x1.*(-2.5e+1./4.0),-x1.*(th13+(th1.*(t15-t21))./4.0e+2),x1.*(th13-(t8.*th1.*(t15-t21))./8.0),t9.*t16.*th1.*x1.*(-2.5e+1./4.0),t10.*t16.*th1.*x1.*(-2.5e+1./4.0),x1.*(t34+th13),t16.*th1.*x1.*(-1.250195899870667e-1),0.0,(t14.*th1.*th18.*x1)./8.0,t6.*t14.*th1.*th18.*x1.*(-1.0./8.0)];
mt5 = [t7.*t14.*th1.*th18.*x1.*(-1.0./8.0),t27,-x1.*(th14+(t8.*th1.*(t13-t14.*x5))./8.0),t9.*t14.*th1.*th18.*x1.*(-1.0./8.0),t10.*t14.*th1.*th18.*x1.*(-1.0./8.0),t27,t14.*th1.*th18.*x1.*(-2.500391799741333e-3),0.0,t56.*th1.*x1.*9.72815332444791e+28,t6.*t56.*th1.*x1.*-9.72815332444791e+28,t7.*t56.*th1.*x1.*-9.72815332444791e+28,t69];
mt6 = [t8.*t56.*th1.*x1.*-9.72815332444791e+28,x1.*(th14-(t9.*th1.*(t51-t51.*t60))./8.0),t10.*t56.*th1.*x1.*-9.72815332444791e+28,t69,t56.*th1.*x1.*(-1.945935583926076e+27),0.0,t50.*th1.*x1.*3.014428576419013e+28];
mt7 = [t6.*t50.*th1.*x1.*-3.014428576419013e+28,t7.*t50.*th1.*x1.*-3.014428576419013e+28,t75,t8.*t50.*th1.*x1.*-3.014428576419013e+28,t9.*t50.*th1.*x1.*-3.014428576419013e+28,x1.*(th14-(t10.*th1.*(t55-t55.*t64))./8.0),t75];
mt8 = [t50.*th1.*x1.*(-6.029801994707232e+26),0.0,t18.*th1.*x1.*(2.5e+1./4.0),t6.*t18.*th1.*x1.*(-2.5e+1./4.0),t7.*t18.*th1.*x1.*(-2.5e+1./4.0),t26,t8.*t18.*th1.*x1.*(-2.5e+1./4.0),t9.*t18.*th1.*x1.*(-2.5e+1./4.0),t10.*t18.*th1.*x1.*(-2.5e+1./4.0),t26,t18.*th1.*x1.*(-1.250195899870667e-1),0.0,t48.*th1.*x1.*6.4101058386717e+28];
mt9 = [t6.*t48.*th1.*x1.*-6.4101058386717e+28,t7.*t48.*th1.*x1.*-6.4101058386717e+28,t70,t8.*t48.*th1.*x1.*-6.4101058386717e+28,t9.*t48.*th1.*x1.*-6.4101058386717e+28,t10.*t48.*th1.*x1.*-6.4101058386717e+28,t70];
mt10 = [t48.*th1.*x1.*(-1.282222085959101e+27),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10],10,10);
