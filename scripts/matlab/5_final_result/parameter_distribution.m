clear ; clc ; close all

experiments = ["all", "bbif", "bbre", "binf", "bthe", "ecol", "lsym", "paci"];

distribution = table();
final_parameters = table();

for i = 1:length(experiments)
	exp_name = experiments(i);
	cd(exp_name)
	new_rows = readtable("distribution.csv");
    distribution = [distribution ; new_rows];
	cd ..
end

for i = 1:length(experiments)
    exp_name = experiments(i);
    cd(exp_name)
    params = readtable("parameters.csv");
    parameter = params.parameter;
    exp = cell(22,1);
    exp(:) = cellstr(exp_name);
    value = params.value;
    new_rows = table(exp, parameter, value);
    final_parameters = [final_parameters ; new_rows];
    cd ..
end

writetable(distribution, "parameter_distribution.csv")
writetable(final_parameters, "parameters.csv")
