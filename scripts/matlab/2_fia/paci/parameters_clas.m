clear ; clc ; close all

% Load parameters
params = readtable("/media/microlab/hdd/taproschle/model/1_first_fit/paci/parameters.csv");

% Separate them in core and shell parameters
core = [1:2, 7:8, 11:14];
shell = [3:6, 9:10, 15:22];

params_core = params(ismember(params.index, core), :);
params_shell = params(ismember(params.index, shell), :);

writetable(params_core, "core_parameters.csv", 'Delimiter', ',')
writetable(params_shell, "shell_parameters.csv", 'Delimiter', ',')