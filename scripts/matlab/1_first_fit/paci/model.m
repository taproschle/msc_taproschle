function dsdt = model(~, s, params)
    dsdt = zeros(10,1);

    % Variables
    X_Com       = s(1);
    S_LNT       = s(2);
    S_2FL       = s(3);
    S_3SL       = s(4);
    S_Lactose   = s(5);
    S_Galactose = s(6);
    S_Glucose   = s(7);
    S_Neu5Ac    = s(8);
    S_Fucose    = s(9);
    S_LNB       = s(10);

    % Parameters
    mu_max = params(1);
    kd     = params(2);

    Y_XCom_SLNT       = params(3);
    Y_XCom_S2FL       = params(4);
    Y_XCom_S3SL       = params(5);
    Y_XCom_SLactose   = params(6);
    Y_XCom_SGalactose = params(7);
    Y_XCom_SGlucose   = params(8);
    Y_XCom_SNeu5Ac    = params(9);
    Y_XCom_SFucose    = params(10);

    k1 = params(11);
    k2 = params(12);
    k3 = params(13);
    k4 = params(14);
    
    K_SLNT       = params(15);
    K_S2FL       = params(16);
    K_S3SL       = params(17);
    K_SLactose   = params(18);
    K_SGalactose = params(19);
    K_SGlucose   = params(20);
    K_SNeu5Ac    = params(21);
    K_SFucose    = params(22);
    
    % Algebraic equations
    mu = mu_max*( ...
        S_LNT/(K_SLNT + S_LNT) + ...
        S_2FL/(K_S2FL + S_2FL) + ...
        S_3SL/(K_S3SL + S_3SL) + ...
        S_Lactose/(K_SLactose + S_Lactose) + ...
        S_Galactose/(K_SGalactose + S_Galactose) + ...
        S_Glucose/(K_SGlucose + S_Glucose) + ...
        S_Neu5Ac/(K_SNeu5Ac + S_Neu5Ac) + ...
        S_Fucose/(K_SFucose + S_Fucose))*(1/8);

    % Equations
    dsdt(1)  =  X_Com*(mu - kd);
    dsdt(2)  = -X_Com*(k1*S_LNT + mu/Y_XCom_SLNT);
    dsdt(3)  = -X_Com*(k2*S_2FL + mu/Y_XCom_S2FL);
    dsdt(4)  = -X_Com*(k3*S_3SL + mu/Y_XCom_S3SL);
    dsdt(5)  =  X_Com*(k1*S_LNT + k2*S_2FL + k3*S_3SL - k4*S_Lactose- mu/Y_XCom_SLactose);
    dsdt(6)  =  X_Com*(k4*S_Galactose - mu/Y_XCom_SGalactose);
    dsdt(7)  =  X_Com*(k4*S_Glucose - mu/Y_XCom_SGlucose);
    dsdt(8)  =  X_Com*(k3*S_3SL - mu/Y_XCom_SNeu5Ac);
    dsdt(9)  =  X_Com*(k2*S_2FL - mu/Y_XCom_SFucose);
    dsdt(10) =  X_Com*k1*S_LNT;

end
