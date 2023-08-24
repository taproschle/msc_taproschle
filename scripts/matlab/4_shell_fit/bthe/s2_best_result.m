% Picks the best results in terms of objective function and SSE
clear ; clc ; close all

% Load data & filter by experiment
data = readtable("/media/microlab/hdd/taproschle/model/reactor_data.csv");
data = data(data.experiment == "bthe", :);
data = data(strcmp(data.unit, 'g/L') | strcmp(data.unit, 'od600'), :);
exclude = {'acetate', 'succinate', 'lactate', 'butirate'};
data = data(~ismember(data.measurement, exclude), :);
pfix = readtable("/media/microlab/hdd/taproschle/model/2_fia/bthe/core_parameters.csv");

% Integrator setup
tspan   = unique(data.time);
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative', 1:10);

% Initial values
y0 = [
    data(data.time == 0 & data.measurement == "biomass", :).value ;
    data(data.time == 0 & data.measurement == "lnt", :).value ;
    data(data.time == 0 & data.measurement == "2fl", :).value ;
    data(data.time == 0 & data.measurement == "3sl", :).value ;
    data(data.time == 0 & data.measurement == "lactose", :).value ;
    data(data.time == 0 & data.measurement == "galactose", :).value ;
    data(data.time == 0 & data.measurement == "glucose", :).value ;
    data(data.time == 0 & data.measurement == "neuac", :).value ;
    data(data.time == 0 & data.measurement == "fucose", :).value ;
    0 ;
    ];

exp_names = ["od600", "lnt", "2fl", "3sl", "lactose", "galactose", ...
    "glucose", "neuac", "fucose"];

fits = zeros(10, 3);

for i = 1:10
    % Load result
    load("results/ess_report_" + num2str(i) + ".mat", "Results")
    k = Results.xbest';
    fun = @(t,y) model(t,y,k,pfix);
    
    % Solve the system
    [tpred, ypred] = ode15s(fun, tspan, y0, options);
    
    % Get SSE of fit
    sse = 0;

    for j = 1:length(exp_names)
        met = exp_names(j);
        met_time = data(data.measurement == met, :).time;
        met_value = data(data.measurement == met, :).value;
        id = zeros(length(met_time), 1);

        for k = 1:length(met_time)
            id(k) = find(tpred == met_time(k));
        end

        met_ypred = ypred(id, j);
        res = met_value - met_ypred;
        sse_j = sum(res.^2);
        sse = sse + sse_j;
    end

    fits(i, :) = [i sse Results.fbest];
end

fits = sortrows(fits, 2);
disp(fits)

best_result_file = "results/ess_report_" + num2str(fits(1,1)) + ".mat";
copyfile(best_result_file, "best_result.mat")