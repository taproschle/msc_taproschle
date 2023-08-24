function dfdxSym = dfdx_model(t,in2,in3)
%DFDX_MODEL
%    dfdxSym = DFDX_MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 22:55:02

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th4 = in3(4,:);
th5 = in3(5,:);
th6 = in3(6,:);
th7 = in3(7,:);
th8 = in3(8,:);
th10 = in3(9,:);
th11 = in3(10,:);
th12 = in3(11,:);
th13 = in3(12,:);
th14 = in3(13,:);
th16 = in3(14,:);
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
t5 = th16+x3;
t6 = th18+x5;
t7 = 1.0./th3;
t8 = 1.0./th4;
t9 = 1.0./th5;
t10 = 1.0./th6;
t11 = 1.0./th7;
t12 = 1.0./th8;
t13 = 1.0./th10;
t14 = x4+5.0e+1;
t28 = x2.*3.5184372088832e+13;
t29 = x7.*1.40737488355328e+14;
t30 = x8.*7.0368744177664e+13;
t31 = x6.*5.62949953421312e+14;
t33 = x9.*9.007199254740992e+15;
t36 = x2+4.99999999947016e+1;
t37 = x7+4.50453941032787e+1;
t38 = x8+9.1442624205459;
t41 = x6+1.43065139909203e+1;
t47 = x9+1.25906487716113e-1;
t15 = 1.0./t5;
t17 = 1.0./t6;
t19 = 1.0./t14;
t35 = t28+1.759218604255179e+15;
t39 = t29+6.339575628071347e+15;
t40 = t30+6.43470262964821e+14;
t42 = t31+8.053851384809931e+15;
t43 = 1.0./t36;
t45 = t33+1.134064822323629e+15;
t50 = 1.0./t37;
t51 = 1.0./t38;
t54 = 1.0./t41;
t56 = 1.0./t47;
t16 = t15.^2;
t18 = t17.^2;
t20 = t19.^2;
t21 = t15.*x3;
t23 = t17.*x5;
t24 = -t15;
t25 = t19.*x4;
t27 = -t19;
t44 = t43.^2;
t46 = 1.0./t35.^2;
t48 = 1.0./t39.^2;
t49 = 1.0./t40.^2;
t52 = t43.*x2;
t55 = 1.0./t42.^2;
t57 = t50.*x7;
t58 = -t43;
t59 = 1.0./t45.^2;
t60 = t51.*x8;
t61 = t54.*x6;
t62 = t56.*x9;
t22 = t16.*x3;
t26 = t20.*x4;
t53 = t44.*x2;
t64 = t21+t23+t25+t52+t57+t60+t61+t62;
t32 = t22+t24;
t34 = t26+t27;
t63 = t53+t58;
mt1 = [-th2+(t64.*th1)./8.0,-t2-(t7.*t64.*th1)./8.0,-t3-(t8.*t64.*th1)./8.0,-t4-(t9.*t64.*th1)./8.0,t2+t3+t4-th14.*x5-(t10.*t64.*th1)./8.0,th14.*x6-(t11.*t64.*th1)./8.0,th14.*x7-(t12.*t64.*th1)./8.0,t4-(t64.*th1)./4.0e+2,t3-(t13.*t64.*th1)./8.0,t2,t46.*th1.*x1.*7.737125244713738e+27,-x1.*(th11+(t7.*th1.*(t43-t53))./8.0)];
mt2 = [t8.*t46.*th1.*x1.*-7.737125244713738e+27,t9.*t46.*th1.*x1.*-7.737125244713738e+27,x1.*(th11-(t10.*th1.*(t43-t53))./8.0),t11.*t46.*th1.*x1.*-7.737125244713738e+27,t12.*t46.*th1.*x1.*-7.737125244713738e+27,t46.*th1.*x1.*(-1.547425048942748e+26)];
mt3 = [t13.*t46.*th1.*x1.*-7.737125244713738e+27,th11.*x1,(t16.*th1.*th16.*x1)./8.0,t7.*t16.*th1.*th16.*x1.*(-1.0./8.0),-x1.*(th12+(t8.*th1.*(t15-t22))./8.0),t9.*t16.*th1.*th16.*x1.*(-1.0./8.0),x1.*(th12-(t10.*th1.*(t15-t22))./8.0),t11.*t16.*th1.*th16.*x1.*(-1.0./8.0),t12.*t16.*th1.*th16.*x1.*(-1.0./8.0),t16.*th1.*th16.*x1.*(-1.0./4.0e+2),x1.*(th12-(t13.*th1.*(t15-t22))./8.0),0.0,t20.*th1.*x1.*(2.5e+1./4.0),t7.*t20.*th1.*x1.*(-2.5e+1./4.0),t8.*t20.*th1.*x1.*(-2.5e+1./4.0)];
mt4 = [-x1.*(th13+(t9.*th1.*(t19-t26))./8.0),x1.*(th13-(t10.*th1.*(t19-t26))./8.0),t11.*t20.*th1.*x1.*(-2.5e+1./4.0),t12.*t20.*th1.*x1.*(-2.5e+1./4.0),x1.*(th13-(th1.*(t19-t26))./4.0e+2),t13.*t20.*th1.*x1.*(-2.5e+1./4.0),0.0,(t18.*th1.*th18.*x1)./8.0,t7.*t18.*th1.*th18.*x1.*(-1.0./8.0),t8.*t18.*th1.*th18.*x1.*(-1.0./8.0),t9.*t18.*th1.*th18.*x1.*(-1.0./8.0),-x1.*(th14+(t10.*th1.*(t17-t18.*x5))./8.0),t11.*t18.*th1.*th18.*x1.*(-1.0./8.0),t12.*t18.*th1.*th18.*x1.*(-1.0./8.0),t18.*th1.*th18.*x1.*(-1.0./4.0e+2),t13.*t18.*th1.*th18.*x1.*(-1.0./8.0),0.0];
mt5 = [t55.*th1.*x1.*5.66739407742615e+29,t7.*t55.*th1.*x1.*-5.66739407742615e+29,t8.*t55.*th1.*x1.*-5.66739407742615e+29,t9.*t55.*th1.*x1.*-5.66739407742615e+29,t10.*t55.*th1.*x1.*-5.66739407742615e+29,x1.*(th14-(t11.*th1.*(t54-t54.*t61))./8.0)];
mt6 = [t12.*t55.*th1.*x1.*-5.66739407742615e+29,t55.*th1.*x1.*(-1.13347881548523e+28),t13.*t55.*th1.*x1.*-5.66739407742615e+29,0.0,t48.*th1.*x1.*1.115269938916765e+29,t7.*t48.*th1.*x1.*-1.115269938916765e+29];
mt7 = [t8.*t48.*th1.*x1.*-1.115269938916765e+29,t9.*t48.*th1.*x1.*-1.115269938916765e+29,t10.*t48.*th1.*x1.*-1.115269938916765e+29,t11.*t48.*th1.*x1.*-1.115269938916765e+29,x1.*(th14-(t12.*th1.*(t50-t50.*t57))./8.0),t48.*th1.*x1.*(-2.230539877833531e+27)];
mt8 = [t13.*t48.*th1.*x1.*-1.115269938916765e+29,0.0,t49.*th1.*x1.*5.660024290063209e+27,t7.*t49.*th1.*x1.*-5.660024290063209e+27,t8.*t49.*th1.*x1.*-5.660024290063209e+27,t9.*t49.*th1.*x1.*-5.660024290063209e+27];
mt9 = [t10.*t49.*th1.*x1.*-5.660024290063209e+27,t11.*t49.*th1.*x1.*-5.660024290063209e+27,t12.*t49.*th1.*x1.*-5.660024290063209e+27,t49.*th1.*x1.*(-1.132004858012642e+26),t13.*t49.*th1.*x1.*-5.660024290063209e+27,0.0];
mt10 = [t59.*th1.*x1.*1.276843477807671e+30,t7.*t59.*th1.*x1.*-1.276843477807671e+30,t8.*t59.*th1.*x1.*-1.276843477807671e+30,t9.*t59.*th1.*x1.*-1.276843477807671e+30,t10.*t59.*th1.*x1.*-1.276843477807671e+30];
mt11 = [t11.*t59.*th1.*x1.*-1.276843477807671e+30,t12.*t59.*th1.*x1.*-1.276843477807671e+30,t59.*th1.*x1.*(-2.553686955615342e+28),t13.*t59.*th1.*x1.*-1.276843477807671e+30,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11],10,10);
