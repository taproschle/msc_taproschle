% Plot results of parameter fitting
clear ; clc ; close all

experiments = ["all", "bbif", "bbre", "binf", "bthe", "ecol", "lsym", "paci"];

sse_table = table();

for i = 1:length(experiments)
    exp = experiments(i);
    params = readtable(exp + "/parameters.csv");
    k = params.value;

    data = readtable("/media/microlab/hdd/taproschle/model/reactor_data.csv");
    data = data(data.experiment == experiments(i), :);
    data = data(strcmp(data.unit, 'g/L') | strcmp(data.unit, 'od600'), :);
    exclude = {'acetate', 'succinate', 'lactate', 'butirate'};
    data = data(~ismember(data.measurement, exclude), :);

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
    
    tspan = linspace(0, 24, 1000);
    options = odeset('NonNegative', 1:10);
    [t, y] = ode15s(@(t, y) model(t,y,k), tspan, y0, options);

    sse = sse_calc(data, k, y0);
    add_row = table(exp, sse);
    sse_table = [sse_table ; add_row];

    result_plot(data, t, y, k, exp)

end

writetable(sse_table, "sse.csv", 'Delimiter', ',')

function sse = sse_calc(data, k, y0)
    % Get SSE of fit
    
    % Integrator setup
    tspan   = unique(data.time);
    options = odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative', 1:10);

    [tpred, ypred] = ode15s(@(t, y) model(t, y, k), tspan, y0, options);

    sse = 0;

    exp_names = ["od600", "lnt", "2fl", "3sl", "lactose", "galactose", ...
    "glucose", "neuac", "fucose"];

    for i = 1:length(exp_names)
        met = exp_names(i);
        met_time = data(data.measurement == met, :).time;
        met_value = data(data.measurement == met, :).value;
        id = zeros(length(met_time), 1);

        for j = 1:length(met_time)
            id(j) = find(tpred == met_time(j));
        end

        met_ypred = ypred(id, i);
        res = met_value - met_ypred;
        sse_j = sum(res.^2);
        sse = sse + sse_j;
    end

end

function [] = result_plot(data, t, y, k, name)

X_Com       = y(:,1);
S_LNT       = y(:,2);
S_2FL       = y(:,3);
S_3SL       = y(:,4);
S_Lactose   = y(:,5);
S_Galactose = y(:,6);
S_Glucose   = y(:,7);
S_Neu5Ac    = y(:,8);
S_Fucose    = y(:,9);
S_LNB       = y(:,10);

mu_max       = k(1);
K_SLNT       = k(15);
K_S2FL       = k(16);
K_S3SL       = k(17);
K_SLactose   = k(18);
K_SGalactose = k(19);
K_SGlucose   = k(20);
K_SNeu5Ac    = k(21);
K_SFucose    = k(22);

mu = mu_max.*( ...
   S_LNT./(K_SLNT + S_LNT) + ...
   S_2FL./(K_S2FL + S_2FL) + ...
   S_3SL./(K_S3SL + S_3SL) + ...
   S_Lactose./(K_SLactose + S_Lactose) + ...
   S_Galactose./(K_SGalactose + S_Galactose) + ...
   S_Glucose./(K_SGlucose + S_Glucose) + ...
   S_Neu5Ac./(K_SNeu5Ac + S_Neu5Ac) + ...
   S_Fucose./(K_SFucose + S_Fucose)).*(1/8);

X_Com_exp       = [data(data.measurement == "od600", :).time, data(data.measurement == "od600", :).value];
S_LNT_exp       = [data(data.measurement == "lnt", :).time, data(data.measurement == "lnt", :).value];
S_2FL_exp       = [data(data.measurement == "2fl", :).time, data(data.measurement == "2fl", :).value];
S_3SL_exp       = [data(data.measurement == "3sl", :).time, data(data.measurement == "3sl", :).value];
S_Lactose_exp   = [data(data.measurement == "lactose", :).time, data(data.measurement == "lactose", :).value];
S_Galactose_exp = [data(data.measurement == "galactose", :).time, data(data.measurement == "galactose", :).value];
S_Glucose_exp   = [data(data.measurement == "glucose", :).time, data(data.measurement == "glucose", :).value];
S_Neu5Ac_exp    = [data(data.measurement == "neuac", :).time, data(data.measurement == "neuac", :).value];
S_Fucose_exp    = [data(data.measurement == "fucose", :).time, data(data.measurement == "fucose", :).value];

% Plots
tiledlayout(3,2)

nexttile
plot(t, X_Com, 'Color', '#003f5c')
grid on
hold on
plot(X_Com_exp(:,1), X_Com_exp(:,2), 's', 'MarkerFaceColor', '#003f5c', 'MarkerEdgeColor', '#003f5c')
legend('Consortium Biomass', 'Consortium Biomass exp.', 'Location', 'best')
xlim([0 24])
xlabel("Time (h)")
ylabel("OD600")

nexttile
plot(t, S_LNT, 'Color', '#003f5c')
grid on
hold on
plot(S_LNT_exp(:,1), S_LNT_exp(:,2), 's', 'MarkerFaceColor', '#003f5c', 'MarkerEdgeColor', '#003f5c')
plot(t, S_2FL, 'Color', '#bc5090')
plot(S_2FL_exp(:,1), S_2FL_exp(:,2), 's', 'MarkerFaceColor', '#bc5090', 'MarkerEdgeColor', '#bc5090')
plot(t, S_3SL, 'Color', '#ffa600')
plot(S_3SL_exp(:,1), S_3SL_exp(:,2), 's', 'MarkerFaceColor', '#ffa600', 'MarkerEdgeColor', '#ffa600')
legend("LNT", "LNT exp.", "2'FL", "2'FL exp.", "3'SL", "3'SL exp.", 'Location','best')
xlim([0 24])
xlabel("Time (h)")
ylabel("g/L")

nexttile
plot(t, mu, 'Color', '#003f5c')
grid on
xlim([0 24])
legend('Growth rate (\mu)')
xlabel("Time (h)")
ylabel("h^{-1}")

nexttile
plot(t, S_Lactose, 'Color', '#003f5c')
grid on
hold on
plot(S_Lactose_exp(:,1), S_Lactose_exp(:,2), 's', 'MarkerFaceColor', '#003f5c', 'MarkerEdgeColor', '#003f5c')
plot(t, S_Galactose, 'Color', '#bc5090')
plot(S_Galactose_exp(:,1), S_Galactose_exp(:,2), 's', 'MarkerFaceColor', '#bc5090', 'MarkerEdgeColor', '#bc5090')
plot(t, S_Glucose, 'Color', '#ffa600')
plot(S_Glucose_exp(:,1), S_Glucose_exp(:,2), 's', 'MarkerFaceColor', '#ffa600', 'MarkerEdgeColor', '#ffa600')
legend("Lactose", "Lactose exp.", "Galactose", "Galactose exp.", "Glucose", "Glucose exp.", 'Location', 'best')
xlim([0 24])
xlabel("Time (h)")
ylabel("g/L")

nexttile
plot(t, S_Neu5Ac, 'Color', '#003f5c')
hold on
grid on
plot(S_Neu5Ac_exp(:,1), S_Neu5Ac_exp(:,2), 's', 'MarkerFaceColor', '#003f5c', 'MarkerEdgeColor', '#003f5c')
plot(t, S_Fucose, 'Color', '#7a5195')
plot(S_Fucose_exp(:,1), S_Fucose_exp(:,2), 's', 'MarkerFaceColor', '#7a5195', 'MarkerEdgeColor', '#7a5195')
legend("Neu5Ac", 'Neu5Ac exp.', "Fucose", 'Fucose exp.', 'Location', 'best')
xlim([0 24])
xlabel("Time (h)")
ylabel("g/L")

name = char(name);
name = [upper(name(1)), lower(name(2:end))];


titleObj = sgtitle(string(name) + " model");
titleObj.FontSize = 16;
titleObj.FontWeight = 'bold';

width = 16;
height = 9;
set(gcf, 'PaperPosition', [0 0 width height])

filename = "plots/" + name + ".png";
print(filename, '-dpng', '-r300')

end