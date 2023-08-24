function f = obj_fun(p)

data = readtable("/media/microlab/hdd/taproschle/model/reactor_data.csv");
data = data(data.experiment == "paci", :);
data = data(strcmp(data.unit, 'g/L') | strcmp(data.unit, 'od600'), :);
exclude = {'acetate', 'succinate', 'lactate', 'butirate'};
data = data(~ismember(data.measurement, exclude), :);

pfix = readtable("/media/microlab/hdd/taproschle/model/3_core_fit/paci/core_parameters.csv");

% Integrator setup
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative', 1:10, ...
    'Events', @timeEvent);

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

od_sta = mean(data(data.measurement == "biomass", :).value(end-2:end));
od_ft = data(data.measurement == "biomass", :).time(end);

time = transpose((od_ft + 1):24);
experiment(1:length(time), 1) = unique(data.experiment);
measurement(1:length(time), 1) = {'biomass'};
unit(1:length(time), 1) = {'od600'};
value(1:length(time), 1) = od_sta;

table_to_add = table(experiment, time, measurement, value, unit);

data = [data ; table_to_add];

tspan   = unique(data.time);

mets = unique(data.measurement);
mets_id = [3 4 1 9 6 7 5 2 8];
weights = [2 2 5 1 1 1 2 2 1];
% weights order
% 2fl 3sl biomass fucose galactose glucose lactose lnt neuac

fun = @(t,y) model(t,y,p,pfix);
tic;
[tpred, ypred] = ode15s(fun, tspan, y0, options);

n_vars = length(y0);
n_times = length(unique(data.time));

% Filter if the integrator found a complete solution
if all(size(ypred) == [n_times n_vars])
    % Init objective function value
    f = 0;
    % Look through all the mets of experimental data
    for i = 1:length(mets)
        % Extract sample times of the mets and values
        met_time = data(data.measurement == string(mets{i}), :).time;
        exp_data = data(data.measurement == string(mets{i}), :).value;
        % Look for the time index for the predicted values
        id = zeros(length(met_time), 1);
        for j = 1:length(met_time)
            id(j) = find(tpred == met_time(j));
        end
        % Extract data
        pred_data = ypred(id, mets_id(i));
        % Add value to the final objective function
        %fi = weights(i)*sum((exp_data - pred_data).^2);
        fi = weights(i)*sum(met_time.*(exp_data - pred_data).^2);
        %fi = weights(i)*sum(abs(exp_data - pred_data));
        %fi = weights(i)*sum(met_time.*abs(exp_data - pred_data));
        f = f + fi;
    end

else
    f = Inf;
end

end
