% Randomly draw new theta's from a uniform distribution

% Define lower and upper bounds for the analysis
thetaLow    = 0.1*thetaNom; 
thetaHigh   = 10*thetaNom;

% Generate the random values for Monte Carlo procedure
THETAReal   = zeros(NExp,nth);

for i=1:nth
    THETAReal(:,i) = random('unif',thetaLow(i),thetaHigh(i),NExp,1);
end

% First draw is nominal value
THETAReal(1,:)  = thetaNom';

% Model Integration and Monte Carlo Procedure

% Define time as a linear space
time    = linspace(0,Tf,Nt);

% Here storage elements are designed for results from integrations during Monte Carlo procedure

dxdth   = cell(NExp,1); %Integrated model jacobian respecting parameters for each iteration 
Xstate  = cell(NExp,1); %Integrated model state's values for each iteration
dxdthRel= cell(NExp,1); %Integrated model relative jacobian respecting parameters for each iteration 

% Additional integration options
options=odeset('RelTol',1e-6,'AbsTol',1e-8);

% Auxiliar functions used in all iterations
X0Model = eval(['@x0_' ModelName]);      %Initial conditions generator for model states
dxdth0  = eval(['@dICdth_',ModelName]);  %Initial conditions generator for matrix dx_th/dt (I.C. for Eq.3 in Stigter & Molenaar (2015))

parfor k=1:NExp
    
    % This iteration's nominal parameter values
    theta=THETAReal(k,:)';
    
    % Generate initial conditions for integration
    x0=X0Model(theta);
    dxdthIC=dxdth0(x0,theta);
        
    % Numerical integration to obtain model dynamics (x_theta in Eq.3 in Stigter & Molenaar (2015))
    [T,Xst]=ode15s(@(t,x) meta(ModelName,t,x,theta,[nx nth]),time,...
       [x0; dxdthIC(:)],options);
   
    % Save columns from 1 to nx corresponding to integrated model states given the inputs
    Xstate{k}=Xst(:,1:nx);
    % Save columns posterior to column nx corresponding to integrated derivative of state vector respect to parameters.
    dxdth{k}=Xst(:,(nx+1):end); 
    
    fprintf('\n Simulation %s done!',num2str(k))
end
