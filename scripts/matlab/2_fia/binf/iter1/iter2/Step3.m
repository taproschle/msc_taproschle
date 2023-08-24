% Initial Procedures for SVD analysis 

% In first place we define those model states acting as sensors (i.e, are measured)
outputsSym = [x1 x2 x3 x4 x5 x6 x7 x8 x9]';       %Substrate

Yobs    = outputsSym;
ny      = length(Yobs);

% Here, gradients are calculated symbolically and automatically stored as
% separated Matlab scripts for later usage. In this case, these are related
% to the output vector.
dhdxSym = simplify(jacobian(Yobs,statesSym));
dhdthSym= simplify(jacobian(Yobs,thetaSym));

h       = matlabFunction(Yobs,'vars',{t,statesSym,thetaSym});
dhdx    = matlabFunction(dhdxSym,'vars',{t,statesSym,thetaSym});
dhdth   = matlabFunction(dhdthSym,'vars',{t,statesSym,thetaSym});

% Here stored elements are designed for results related to parametric
% output sensitivities
dydth   = cell(NExp,1);  %Partial derivates of the output vector with respect theta.
dydthRel= cell(NExp,1);  %Relative partial derivates of the output vector with respect theta.
YModel  = cell(NExp,1);  %Integrated measured model states 

% Construction of the ROSM for each iteration of Monte Carlo procedure
for k=1:NExp
    
    % Here we construct the ROSM (Stigter & Molenaar, 2015):
    
    dydth{k}    = zeros(ny*Nt,nth);
    dydthRel{k} = zeros(ny*Nt,nth);
    Ym          = zeros(Nt,ny);

    for i=1:Nt
        Ym(i,:) = h(time(i),Xstate{k}(i,1:nx)',THETAReal(k,:)')';
    end
    YModel{k}   = Ym;

    for i=1:Nt
           
        dhdxM   = dhdx(time(i),Xstate{k}(i,1:nx)',thetaNom);
        dhdthM  = dhdth(time(i),Xstate{k}(i,1:nx)',thetaNom);
        dydthMatrix = dhdxM*reshape(dxdth{k}(i,:)',nx,nth)+dhdthM;
        
        % Observation Matrix Normalization
        dydthRelMatrix = dydthMatrix.*(ones(ny,1)*thetaNom').*  ...
            ((1./(Ym(i,:)')*ones(1,nth)));
        
        % Iteration specific definitions and storage variable actualizations 
        Index               = ((i-1)*ny+1):(i*ny);
        dydth{k}(Index,:)   = dydthMatrix;
        dydthRel{k}(Index,:)= dydthRelMatrix;
    end
end

% Perform SVD and store results

% Define storage variables
svdResults          = cell(NExp,3);
svdResultsAbsSens   = cell(NExp,3);
SingularValues      = zeros(NExp,nth);
NormRelSen          = cell(NExp,1);
RankingPar          = cell(NExp,1);

% SVD decomposition cycle for each k-th Monte Carlo iteration 
for k=1:NExp
    [U,S,V]             = svd(dydthRel{k});
    svdResults{k,1}     = U;
    svdResults{k,2}     = S;
    svdResults{k,3}     = V;
    SingularValues(k,:) = diag(S)';
    NormRelSen{k}       = sum(dydth{k}.*dydth{k});
end

% Optional if some parameters are to be discarded in the analysis
ParIndex    = 1:nth;
ExcludePar  = [];
ParIndex    = setdiff(ParIndex,ExcludePar);
nthA        = numel(ParIndex);

% Graphical Artwork

% Figure 1 corresponding to the values of each component of the V vector
% associated with the smallest singular value (S) from the SVD
% decomposition of the ROSM in each k-th iteration

figure(1);
hold on

% Find smallest singular values and create stored variables
minS    = min(SingularValues(:,end));
v_vec   = zeros(nth,NExp);

% Stem plot of V vector components
for k = 1:NExp
    V = svdResults{k,3};
    stem(V(:,nth));
    v_vec(:,k) = abs(V(:,nth));
end

% Save statistics of the component values for the NExp experiments
mean_vec = mean(v_vec,2);
std_vect = std(v_vec,0,2);

% Settings for the lables in the x axis
set(gca,'XTick',1:nth);
xtick=cell(nth,1);
for i=1:nth
    xtick{i}=char(thetaSym(ParIndex(i)));
end
set(gca,'XTickLabel',xtick);

% Other settings
textbp(['\sigma_{',num2str(nth),'}=',num2str(minS,'%10.5e\n')],'Fontname','FixedWidth');
title('Last column of V (Nullspace)','FontWeight','bold');
axis([0.7,nth+0.3,-1,1]);
xlabel('Parameter');
ylabel('Component Last Singular Vector');


% Figure 2 corresponding to the mean of the absolut values obtained from
% the Monte Carlo simulations associated to each component of vector V

figure(2)

% Stem plot for the mean values obtained from the Monte Carlo procedure
% (mean absolute values of components associated to V singular vector for
% each iteration)
stem(mean_vec);

% Settings for the lables in the x axis
set(gca,'XTick',1:nth);
xtick=cell(nth,1);
for i=1:nth
    xtick{i}=char(thetaSym(ParIndex(i)));
end
set(gca,'XTickLabel',xtick);

%  Other settings
textbp(['\sigma_{',num2str(nth),'}=',num2str(minS,'%10.5e\n')],'Fontname','FixedWidth');
title('Last column of V (Nullspace)','FontWeight','bold');
axis([0.7,nth+0.3,-1,1]);
xlabel('Parameter');
ylabel('Component Last Singular Vector');
lblc = sprintfc('%.2f +/- %.2f',[mean_vec,std_vect]);
text(1:nth, mean_vec+0.05, lblc, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')

% Figure 3 corresponding to the plotting of the singular values obtained in
% each iteration

% Save min and maximum singular values obtained from the analysis
minS=min(SingularValues(:,end));
maxS=max(SingularValues(:,1));

figure(3);
plot(log10(SingularValues),'.','MarkerSize',18);
xlabel('Virtual Experiment Number');
ylabel('^{10}Log(SingularValues)');
title('10Log(Singular Values)','FontWeight','bold');
axis([0 NExp+1 floor(log10(minS)) ceil(log10(maxS))]);

