% Runs parameter fitting

x_B = [ 0.1     2;   % mu_max
        1e-5    0.1; % kd
        1e-5    75;  % Y_XCom_SLNT
        1e-5    75;  % Y_XCom_S2FL
        1e-5    75;  % Y_XCom_S3SL
        1e-5    75;  % Y_XCom_SLactose
        1e-5    75;  % Y_XCom_SGalactose
        1e-5    75;  % Y_XCom_SGlucose
        1e-5    75;  % Y_XCom_SNeu5Ac
        1e-5    75;  % Y_XCom_SFucose
        1e-5    75;  % k1
        1e-5    75;  % k2
        1e-5    75;  % k3
        1e-5    75;  % k4
        1e-5    75;  % K_SLNT
        1e-5    75;  % K_S2FL
        1e-5    75;  % K_S3SL
        1e-5    75;  % K_SLactose
        1e-5    75;  % K_SGalactose
        1e-5    75;  % K_SGlucose
        1e-5    75;  % K_SNeu5Ac
        1e-5    75]; % K_SFucose

pfix = readtable("/media/microlab/hdd/taproschle/model/3_core_fit/ecol/core_parameters.csv");
x_B(pfix.index, :) = [];

% MEIGO Settings
problem.f   = 'obj_fun';
problem.x_L = x_B(:,1);
problem.x_U = x_B(:,2);

opts.maxeval = 10000;
opts.maxtime = 500;
opts.iterprint =  1;

opts.ndiverse = 2000;
opts.local.solver = 'fmincon';
opts.local.finish = 'fmincon';

opts.n_threads = 8;
opts.n_iter = 1;
opts.is_parallel = true;
opts.maxtime_per_iteration = 10;
opts.maxtime_cess = 300;
opts.log_var = [];

id = 'MATLAB:ode15s:IntegrationTolNotMet';
warning('off',id)

wd = pwd();

dir = wd + "/results";

if ~exist(dir, 'dir')
    mkdir(dir);
end

for i = 1:10
    % Pick random starting point
    problem.x_0 = x_B(:,1) + rand(size(x_B,1),1) .* (x_B(:,2) - x_B(:,1));
    % Run scatter search
    MEIGO(problem, opts, 'CESS');
    % Save file
    file = wd + "/ess_report.mat";
    file_new = "ess_report_" + num2str(i) + ".mat";
    file_new = fullfile(dir, file_new);
    movefile(file, file_new)
end
