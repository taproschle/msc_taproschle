function dfdthSym = dfdth_model(t,in2,in3)
%DFDTH_MODEL
%    dfdthSym = DFDTH_MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 22:19:47

th1 = in3(1,:);
th3 = in3(3,:);
th4 = in3(4,:);
th5 = in3(5,:);
th6 = in3(6,:);
th7 = in3(7,:);
th8 = in3(8,:);
th10 = in3(9,:);
th18 = in3(14,:);
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
t12 = 1.0./th10;
t16 = x4+4.99989432222649e+1;
t17 = x7+2.09075062954877e+1;
t18 = x8+3.31058128076456e+1;
t19 = x3+4.96288289046932e+1;
t20 = x6+4.952428510331e+1;
t21 = x2+1.7925615324441e+1;
t22 = x9+6.18782139387954;
t13 = 1.0./t2;
t23 = 1.0./t16;
t24 = 1.0./t21;
t25 = 1.0./t17;
t26 = 1.0./t18;
t27 = 1.0./t19;
t28 = 1.0./t20;
t29 = 1.0./t22;
t14 = t13.^2;
t15 = t13.*x5;
t30 = t23.*x4;
t31 = t24.*x2;
t32 = t25.*x7;
t33 = t26.*x8;
t34 = t27.*x3;
t35 = t28.*x6;
t36 = t29.*x9;
t37 = t15+t30+t31+t32+t33+t34+t35+t36;
et1 = t15.*2.50000001290012e-3+t30.*2.50000001290012e-3+t31.*2.50000001290012e-3+t32.*2.50000001290012e-3+t33.*2.50000001290012e-3+t34.*2.50000001290012e-3;
et2 = t35.*2.50000001290012e-3+t36.*2.50000001290012e-3;
mt1 = [x1.*(t15./8.0+t30./8.0+t31./8.0+t32./8.0+t33./8.0+t34./8.0+t35./8.0+t36./8.0),t6.*t37.*x1.*(-1.0./8.0),t7.*t37.*x1.*(-1.0./8.0),t8.*t37.*x1.*(-1.0./8.0),t9.*t37.*x1.*(-1.0./8.0),t10.*t37.*x1.*(-1.0./8.0),t11.*t37.*x1.*(-1.0./8.0),-x1.*(et1+et2),t12.*t37.*x1.*(-1.0./8.0),0.0,-x1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t6.^2.*t37.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t7.^2.*t37.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t8.^2.*t37.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
mt2 = [(t9.^2.*t37.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t10.^2.*t37.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t11.^2.*t37.*th1.*x1)./8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(t12.^2.*t37.*th1.*x1)./8.0,0.0,0.0,-t3,0.0,0.0,t3,0.0,0.0,0.0,0.0,t3,0.0,0.0,-t4,0.0,t4,0.0,0.0,0.0,t4,0.0,0.0,0.0,0.0,-t5,t5,0.0,0.0,t5,0.0,0.0,0.0,0.0,0.0,0.0,-x1.*x5,x1.*x6,x1.*x7,0.0,0.0,0.0,t14.*th1.*x1.*x5.*(-1.0./8.0),(t6.*t14.*th1.*x1.*x5)./8.0,(t7.*t14.*th1.*x1.*x5)./8.0,(t8.*t14.*th1.*x1.*x5)./8.0,(t9.*t14.*th1.*x1.*x5)./8.0,(t10.*t14.*th1.*x1.*x5)./8.0];
mt3 = [(t11.*t14.*th1.*x1.*x5)./8.0,t14.*th1.*x1.*x5.*2.50000001290012e-3,(t12.*t14.*th1.*x1.*x5)./8.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3],10,14);
