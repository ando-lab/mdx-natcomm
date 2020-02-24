%%
% load reference data and lattice model

load model/lattice_dynamics_model.mat M N
load fit/reference_halos.mat hklGrid I sigma
load fit/reference_calc.mat K1 Wk

%% 
% initialize variables to store refinement info (4 stages)
history = cell(4,1);
xfit = cell(4,1);
fitinfo = cell(4,1);

%% Stage 1
% fit overall Gaussian spring constant

N.springType = 'gaussian';
optfun = @(x) modelFitFunction(I,sigma,Wk,K1,N,M,x);
x0 = 0.5;
xmin = 0;
xmax = Inf;
fitoptions = {'Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',1000,'UseParallel',true,'Display','iter'};

% run the optimization (wrapper for lsqnonlin)
[history{1},xfit{1},fitinfo{1}] = runlsqnonlin(optfun,x0,xmin,xmax,fitoptions{:});

% assign spring constants
for j=1:numel(history{1})
    history{1}(j).springType = N.springType;
    history{1}(j).springConstants = history{1}(j).x;
end

%% Stage 2
% fit Gaussian spring constant for each interface

[uniqueInterfaces,~,ix] = unique(abs(N.Edges.interface));
gammaWrap = @(g0) g0(ix);
N.springType = 'gaussian';

optfun = @(x) modelFitFunction(I,sigma,Wk,K1,N,M,gammaWrap(x));

x0 = repmat(xfit{1},length(uniqueInterfaces),1);
xmin = x0*0;
xmax = x0*Inf;

fitoptions = {'Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',1000,'UseParallel',true,'Display','iter'};


% run the optimization (wrapper for lsqnonlin)
[history{2},xfit{2},fitinfo{2}] = runlsqnonlin(optfun,x0,xmin,xmax,fitoptions{:});

% assign spring constants
for j=1:numel(history{2})
    history{2}(j).springType = N.springType;
    history{2}(j).springConstants = gammaWrap(history{2}(j).x);
end

%% Stage 3
% fit hybrid spring constants for each interface

[uniqueInterfaces,~,ix] = unique(abs(N.Edges.interface));
gammaWrap = @(g0) g0(ix,:);
N.springType = 'hybrid';

optfun = @(x) modelFitFunction(I,sigma,Wk,K1,N,M,gammaWrap(x));

x0 = repmat(xfit{2},1,2);
xmin = x0*0;
xmax = x0*Inf;

fitoptions = {'Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',1000,'UseParallel',true,'Display','iter'};

% run the optimization (wrapper for lsqnonlin)
[history{3},xfit{3},fitinfo{3}] = runlsqnonlin(optfun,x0,xmin,xmax,fitoptions{:});

% assign spring constants
for j=1:numel(history{3})
    history{3}(j).springType = N.springType;
    history{3}(j).springConstants = gammaWrap(history{3}(j).x);
end
%% Stage 4
% fit individual hybrid springs
T = N.Edges;
ifswitch = T.interface < 1;
SpringInfo = T;
SpringInfo.a1(ifswitch) = T.a2(ifswitch);
SpringInfo.a2(ifswitch) = T.a1(ifswitch);
SpringInfo.interface(ifswitch) = -T.interface(ifswitch);

[uniqueSprings,ix0,ix] = unique(SpringInfo(:,{'interface','a1','a2'}),'rows');
gfit0 = gammaWrap(xfit{3});
x0 = gfit0(ix0,:);

gammaWrap = @(g0) g0(ix,:);

N.springType = 'hybrid';

optfun = @(x) modelFitFunction(I,sigma,Wk,K1,N,M,gammaWrap(x));

xmin = x0*0;
xmax = x0*Inf;

fitoptions = {'Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',3000,'UseParallel',true,'Display','iter'};

% run the optimization (wrapper for lsqnonlin)
[history{4},xfit{4},fitinfo{4}] = runlsqnonlin(optfun,x0,xmin,xmax,fitoptions{:});

% assign spring constants
for j=1:numel(history{4})
    history{4}(j).springType = N.springType;
    history{4}(j).springConstants = gammaWrap(history{4}(j).x);
end

%% save result

N.springConstants = gammaWrap(xfit{4});
K = N.Hessian();
L = nm.LatticeModel(M,K,K1.N);

save fit/fit_lattice_dynamics_model_to_reference.mat L N history xfit fitinfo

%% CONSOLE OUTPUT
% (Feb 24, 2020)
%
%
%
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 10 workers.
% 
%                                          Norm of      First-order 
%  Iteration  Func-count     f(x)          step          optimality
%      0          2     2.84986e+06                      1.78e+06
%      1          4     2.48747e+06       0.119808         0.0844      
%      2          6     2.48747e+06    4.31015e-09         0.0844      
% 
% Local minimum possible.
% 
% lsqnonlin stopped because the size of the current step is less than
% the default value of the step size tolerance.
% 
% <stopping criteria details>
% 
% 
%                                          Norm of      First-order 
%  Iteration  Func-count     f(x)          step          optimality
%      0          7     2.48747e+06                      5.38e+05
%      1         14     2.14968e+06       0.616504       2.92e+05      
%      2         21     2.10915e+06        0.33071       4.39e+04      
%      3         28     2.10859e+06      0.0134432       2.06e+03      
%      4         35     2.10858e+06     0.00616004            167      
%      5         42     2.10858e+06    0.000316021           39.8      
% 
% Local minimum possible.
% 
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the default value of the function tolerance.
% 
% <stopping criteria details>
% 
% 
%                                          Norm of      First-order 
%  Iteration  Func-count     f(x)          step          optimality
%      0         13     2.10858e+06                      8.44e+04
%      1         26     1.92289e+06       0.994418       1.61e+05      
%      2         39     1.91581e+06      0.0298967       3.42e+04      
%      3         52      1.8829e+06       0.427908       2.83e+04      
%      4         65     1.87034e+06       0.279809       1.42e+04      
%      5         78     1.86569e+06       0.198595       5.08e+03      
%      6         91     1.86402e+06       0.148272       2.29e+03      
%      7        104     1.86339e+06       0.113963       1.85e+03      
%      8        117     1.86316e+06      0.0860797       1.25e+03      
%      9        130     1.86308e+06      0.0628586            785      
%     10        143     1.86305e+06      0.0440499            476      
%     11        156     1.86304e+06      0.0289832            277      
%     12        169     1.86304e+06      0.0172399            150      
%     13        182     1.86304e+06     0.00892309           73.5      
% 
% Local minimum possible.
% 
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the default value of the function tolerance.
% 
% <stopping criteria details>
% 
% 
%                                          Norm of      First-order 
%  Iteration  Func-count     f(x)          step          optimality
%      0        101     1.86304e+06                      3.13e+04
%      1        202      1.8217e+06       0.405815       2.51e+04      
%      2        303     1.57292e+06        7.56624       4.11e+04      
%      3        404     1.55462e+06      0.0922331       1.41e+04      
%      4        505     1.54947e+06       0.113516       7.02e+03      
%      5        606     1.52439e+06        1.75395       6.94e+03      
%      6        707     1.48802e+06        2.04464       2.63e+03      
%      7        808     1.46563e+06        3.90195       1.72e+03      
%      8        909     1.45006e+06        2.63388            955      
%      9       1010     1.44362e+06         1.7466            586      
%     10       1111     1.44135e+06        1.01138       3.91e+03      
%     11       1212     1.43918e+06       0.529983            187      
%     12       1313     1.43845e+06       0.771257       3.74e+03      
%     13       1414     1.43788e+06       0.308375           64.8      
%     14       1515     1.43754e+06       0.758209             72      
%     15       1616     1.43747e+06        0.30791           54.4      
%     16       1717     1.43743e+06       0.196272           47.5      
%     17       1818      1.4374e+06       0.177829           39.3      
%     18       1919     1.43739e+06       0.154146           34.4      
%     19       2020     1.43738e+06       0.124979           29.5      
%     20       2121     1.43738e+06      0.0947799           23.9      
%     21       2222     1.43737e+06      0.0543098           17.6      
%     22       2323     1.43737e+06      0.0290765           11.8      
%
% Local minimum possible.
% 
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the default value of the function tolerance.
% 
% <stopping criteria details>
%