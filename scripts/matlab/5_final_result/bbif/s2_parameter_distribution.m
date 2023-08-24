clear; clc; close all

% Load parameters
core = readtable("../../3_core_fit/bbif/core_parameters.csv");
shell = readtable("../../4_shell_fit/bbif/shell_parameters.csv");

% Create empty table
exp = {'bbif'};
parameter_distribution = table();

% Search for core parameters values
for i = 1:height(core)
    for j = 1:20
        load("../../3_core_fit/bbif/results/ess_report_" + num2str(j) + ".mat", 'Results')
        parameter = core.parameter(i);
        value = Results.xbest(i);
        new_row = table(exp, parameter, value);
        parameter_distribution = [parameter_distribution; new_row];
    end
end

% Search for shell parameters values
for i = 1:height(shell)
    for j = 1:10
        load("../../4_shell_fit/bbif/results/ess_report_" + num2str(j) + ".mat", 'Results')
        parameter = shell.parameter(i);
        value = Results.xbest(i);
        new_row = table(exp, parameter, value);
        parameter_distribution = [parameter_distribution; new_row];
    end
end

writetable(parameter_distribution, 'distribution.csv', 'Delimiter', ',')
