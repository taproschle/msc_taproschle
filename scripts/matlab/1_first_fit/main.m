clear ; clc ; close all

experiments = ["all", "bbif", "bbre", "binf", "bthe", "ecol", "lsym", "paci"];

% This part of the workflow can be executed without supervision
for i = 1:length(experiments)
	exp = experiments(i);
	cd(exp)
	run s1_meigo_cess.m
	cd ..
end

% The following scripts s2 and s3 needs to be visually evaluated
