clear ; clc ; close all

% Load parameters
params = readtable("/media/microlab/hdd/taproschle/model/1_first_fit/all/parameters.csv");

% Separate them in core and shell parameters
core = [1:4, 6:8, 12:14];
shell = [5, 9:11, 15:22];

params_core = params(ismember(params.index, core), :);
params_shell = params(ismember(params.index, shell), :);

writetable(params_core, "core_parameters.csv", 'Delimiter', ',')
writetable(params_shell, "shell_parameters.csv", 'Delimiter', ',')