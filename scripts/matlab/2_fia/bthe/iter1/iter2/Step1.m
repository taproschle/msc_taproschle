% Adapted from https://sourceforge.net/projects/minimal-output-sets/
% by Cristóbal Torrealba 2021

% Initial procedures
clear ; clc ; close all
tic
ModelName = 'model';

% Load parameters from previous step (Initial guess)
params = readtable("/media/microlab/hdd/taproschle/model/1_first_fit/bthe/parameters.csv");
k = params.value;

data = readtable("/media/microlab/hdd/taproschle/model/reactor_data.csv");
data = data(data.experiment == "bthe", :);
data = data(strcmp(data.unit, 'g/L') | strcmp(data.unit, 'od600'), :);
exclude = {'acetate', 'succinate', 'lactate', 'butirate'};
data = data(~ismember(data.measurement, exclude), :);

% Initial conditions for integration
x0 = [
    data(data.time == 0 & data.measurement == "od600", :).value ;
    data(data.time == 0 & data.measurement == "lnt", :).value ;
    data(data.time == 0 & data.measurement == "2fl", :).value ;
    data(data.time == 0 & data.measurement == "3sl", :).value ;
    data(data.time == 0 & data.measurement == "lactose", :).value ;
    data(data.time == 0 & data.measurement == "galactose", :).value ;
    data(data.time == 0 & data.measurement == "glucose", :).value ;
    data(data.time == 0 & data.measurement == "neuac", :).value ;
    data(data.time == 0 & data.measurement == "fucose", :).value ;
    0 ;
    ]';

x0 = x0 + 1e-10;

% Define final integration time (Tf, in hours for this example) and
% number of integration steps (Nt) 
Tf = 24;
Nt = 100;

% Symbolic Toolbox Variable defintions

% Here, we define model states as x{i} for a total of i = 1..10 model states.
for i = 1:10
    syms(sprintf('x%d',i),'real');
end

% Here, we define model parameters as th{i} for a total of i = 1..22 regression parameters.
for i = 1:22
    syms(sprintf('th%d',i),'real');
end

% Additionally, we define variables t symbolic values representing time and model inputs.
syms t real

% Once symbolic parameters and model states are defined, we build arrays to vectorize parameters and model states.
statesSym=[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10]';
thetaSym=[th1  th2  th3  th6  th7  th8  th11 th12 ...
          th13 th14]';

% Fixed parameters
th9 = k(9);
th16 = k(16);
th17 = k(17);
th18 = k(18);
th19 = k(19);
th20 = k(20);
th21 = k(21);
th22 = k(22);

th4 = k(4);
th5 = k(5);
th10 = k(10);
th15 = k(15);

k = k([1:3, 6:8, 11:14]);

% Define the amount of model states and parameters
nx=length(statesSym);
nth=length(thetaSym);

% Define initial conditions vector as a symbolic variable associated to numerical values.
thIC=sym(x0);
x0Sym=thIC';

% Define the amount of simulations performed in the Monte Carlo procedure
NExp=100; 

% Nominal parameter values

thetaNom = k;

% MODEL DEFINITION

mu = th1*(x2/(th15 + x2) +  x3/(th16 + x3) +  x4/(th17 + x4) +  x5/(th18 + x5) +  x6/(th19 + x6) +  x7/(th20 + x7) +  x8/(th21 + x8) +  x9/(th22 + x9))*(1/8);

Xdot =  [
    x1*(mu - th2);
    -x1*(th11*x2 + mu/th3);
    -x1*(th12*x3 + mu/th4);
    -x1*(th13*x4 + mu/th5);
    x1*(th11*x2 + th12*x3 + th13*x4 - th14*x5- mu/th6);
    x1*(th14*x6 - mu/th7);
    x1*(th14*x7 - mu/th8);
    x1*(th13*x4 - mu/th9);
    x1*(th12*x3 - mu/th10);
    x1*th11*x2;
    ];


% END MODEL DEFINITION AND AUXILIAR FUNCTION GENERATION
% (May take a few minutes to process)

% Here, gradients are calculated symbolically and automatically stored as separated Matlab scripts for later usage.
dfdxSym     = simplify(jacobian(Xdot,statesSym));
dfdthSym    = simplify(jacobian(Xdot,thetaSym));
dx0dthSym   = simplify(jacobian(x0Sym,thetaSym));

f = matlabFunction(Xdot,'vars',{t,statesSym,thetaSym},'file',ModelName);
dfdx = matlabFunction(dfdxSym,'vars',{t,statesSym,thetaSym},'file',['dfdx_',ModelName]);
dfdth = matlabFunction(dfdthSym,'vars',{t,statesSym,thetaSym},'file',['dfdth_',ModelName]);
IC = matlabFunction(x0Sym,'vars',{thetaSym},'file',['x0_',ModelName]);
dICdth = matlabFunction(dx0dthSym,'vars',{statesSym,thetaSym},'file',['dICdth_',ModelName]);
