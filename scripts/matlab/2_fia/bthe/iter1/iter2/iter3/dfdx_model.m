function dfdxSym = dfdx_model(t,in2,in3)
%DFDX_MODEL
%    dfdxSym = DFDX_MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    31-May-2023 22:47:18

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
t5 = 1.0./th3;
t6 = 1.0./th7;
t7 = 1.0./th8;
t8 = x3.*3.5184372088832e+13;
t9 = x8.*3.5184372088832e+13;
t10 = x4.*1.40737488355328e+14;
t11 = x6.*7.0368744177664e+13;
t12 = x5.*2.251799813685248e+15;
t13 = x7.*4.503599627370496e+15;
t14 = x9.*2.251799813685248e+15;
t15 = x2.*9.007199254740992e+15;
t16 = x8+3.7995870797099e+1;
t19 = x3+4.99999999906868e+1;
t22 = x6+2.19492689619356e+1;
t23 = x4+4.99999994988811e+1;
t25 = x9+1.86124538346878;
t26 = x7+1.66530861806818;
t27 = x5+3.84036420636605;
t34 = x2+8.16450007953423e-1;
t17 = t8+1.759218604113921e+15;
t18 = t9+1.336860855964317e+15;
t20 = t11+1.544542492469187e+15;
t21 = t10+7.036874347240185e+15;
t24 = t12+8.647731404378567e+15;
t28 = t14+4.191152007717527e+15;
t29 = t13+7.499883271788731e+15;
t30 = 1.0./t19;
t32 = 1.0./t16;
t36 = t15+7.353927903171349e+15;
t37 = 1.0./t22;
t38 = 1.0./t23;
t44 = 1.0./t27;
t49 = 1.0./t25;
t50 = 1.0./t26;
t56 = 1.0./t34;
t31 = t30.^2;
t33 = 1.0./t17.^2;
t35 = 1.0./t18.^2;
t39 = t38.^2;
t40 = 1.0./t20.^2;
t41 = 1.0./t21.^2;
t42 = t30.*x3;
t45 = 1.0./t28.^2;
t46 = t32.*x8;
t47 = 1.0./t29.^2;
t48 = 1.0./t24.^2;
t51 = t38.*x4;
t53 = 1.0./t36.^2;
t54 = -t30;
t55 = t37.*x6;
t57 = t56.^2;
t58 = -t38;
t59 = t49.*x9;
t60 = t50.*x7;
t61 = t44.*x5;
t62 = t56.*x2;
t64 = -t56;
t43 = t31.*x3;
t52 = t39.*x4;
t63 = t57.*x2;
t65 = t33.*th1.*x1.*1.547425049106725e+26;
t70 = t35.*th1.*x1.*1.175915244900105e+26;
t72 = t40.*th1.*x1.*2.717187888608508e+26;
t74 = t45.*th1.*x1.*2.359408827965693e+28;
t75 = t48.*th1.*x1.*4.868239992201709e+28;
t77 = t47.*th1.*x1.*8.444117878610321e+28;
t80 = t53.*th1.*x1.*1.655957348530044e+29;
t83 = t42+t46+t51+t55+t59+t60+t61+t62;
t66 = -t65;
t67 = t43+t54;
t68 = t52+t58;
t69 = t63+t64;
t71 = -t70;
t73 = -t72;
t76 = -t74;
t78 = -t75;
t79 = -t77;
t81 = -t80;
t82 = th1.*(t38-t52).*(-2.50000000046566e-3);
t84 = t83.*th1.*2.50000000046566e-3;
t85 = -t84;
mt1 = [-th2+(t83.*th1)./8.0,-t2-(t5.*t83.*th1)./8.0,-t3-t83.*th1.*4.223924006036696e-2,-t4+t85,t2+t3+t4-t83.*th1.*4.153991908225623e-1-th14.*x5,th14.*x6-(t6.*t83.*th1)./8.0,th14.*x7-(t7.*t83.*th1)./8.0,t4+t85,t3-t83.*th1.*2.593753788511239e-1,t2];
mt2 = [t53.*th1.*x1.*8.279786741107995e+30,-x1.*(th11+(t5.*th1.*(t56-t63))./8.0),t53.*th1.*x1.*(-2.797855198450432e+30),t81,x1.*(th11-th1.*(t56-t63).*4.153991908225623e-1)];
mt3 = [t6.*t53.*th1.*x1.*-8.279786741107995e+30,t7.*t53.*th1.*x1.*-8.279786741107995e+30,t81,t53.*th1.*x1.*(-1.718058258225119e+31),th11.*x1,t33.*th1.*x1.*7.737125244092479e+27];
mt4 = [t5.*t33.*th1.*x1.*-7.737125244092479e+27,-x1.*(th12+th1.*(t30-t43).*4.223924006036696e-2),t66,x1.*(th12-th1.*(t30-t43).*4.153991908225623e-1),t6.*t33.*th1.*x1.*-7.737125244092479e+27,t7.*t33.*th1.*x1.*-7.737125244092479e+27,t66];
mt5 = [x1.*(th12-th1.*(t30-t43).*2.593753788511239e-1),0.0,t41.*th1.*x1.*1.237940026878277e+29,t5.*t41.*th1.*x1.*-1.237940026878277e+29,t41.*th1.*x1.*(-4.183171678051895e+28)];
mt6 = [-x1.*(th13+th1.*(t38-t52).*2.50000000046566e-3),x1.*(th13-th1.*(t38-t52).*4.153991908225623e-1),t6.*t41.*th1.*x1.*-1.237940026878277e+29,t7.*t41.*th1.*x1.*-1.237940026878277e+29,x1.*(t82+th13)];
mt7 = [t41.*th1.*x1.*(-2.568729307732189e+29),0.0,t48.*th1.*x1.*2.434119995647466e+30,t5.*t48.*th1.*x1.*-2.434119995647466e+30];
mt8 = [t48.*th1.*x1.*(-8.225230306551415e+29),t78,-x1.*(th14+th1.*(t44-t44.*t61).*4.153991908225623e-1),t6.*t48.*th1.*x1.*-2.434119995647466e+30,t7.*t48.*th1.*x1.*-2.434119995647466e+30,t78];
mt9 = [t48.*th1.*x1.*(-5.050806368321259e+30),0.0,t40.*th1.*x1.*1.358593944051197e+28,t5.*t40.*th1.*x1.*-1.358593944051197e+28];
mt10 = [t40.*th1.*x1.*(-4.590878059787141e+27),t73,t40.*th1.*x1.*(-4.514870600122405e+28),x1.*(th14-(t6.*th1.*(t37-t37.*t55))./8.0),t7.*t40.*th1.*x1.*-1.358593944051197e+28,t73];
mt11 = [t40.*th1.*x1.*(-2.819086551544974e+28),0.0,t47.*th1.*x1.*4.222058938518743e+30,t5.*t47.*th1.*x1.*-4.222058938518743e+30];
mt12 = [t47.*th1.*x1.*(-1.426692488424891e+30),t79,t47.*th1.*x1.*(-1.403071893332682e+31),t6.*t47.*th1.*x1.*-4.222058938518743e+30,x1.*(th14-(t7.*th1.*(t50-t50.*t60))./8.0),t79];
mt13 = [t47.*th1.*x1.*(-8.760785093680584e+30),0.0,t35.*th1.*x1.*5.879576223405371e+27,t5.*t35.*th1.*x1.*-5.879576223405371e+27];
mt14 = [t35.*th1.*x1.*(-1.986790652429162e+27),t71,t35.*th1.*x1.*(-1.953896964465735e+28),t6.*t35.*th1.*x1.*-5.879576223405371e+27];
mt15 = [t7.*t35.*th1.*x1.*-5.879576223405371e+27,t71,t35.*th1.*x1.*(-1.220013848343863e+28),0.0,t45.*th1.*x1.*1.17970441376311e+30,t5.*t45.*th1.*x1.*-1.17970441376311e+30];
mt16 = [t45.*th1.*x1.*(-3.986385434657159e+29),t76,t45.*th1.*x1.*(-3.920386071096009e+30),t6.*t45.*th1.*x1.*-1.17970441376311e+30];
mt17 = [t7.*t45.*th1.*x1.*-1.17970441376311e+30,t76,t45.*th1.*x1.*(-2.447890234017197e+30),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11,mt12,mt13,mt14,mt15,mt16,mt17],10,10);