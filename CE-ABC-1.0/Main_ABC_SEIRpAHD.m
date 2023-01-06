% -----------------------------------------------------------------
%  Main_ABC_SEIRpAHD.m
% -----------------------------------------------------------------
%  This program uses the CE-ABC algorithm to identify parameters
%  and propagate uncertainties in an SEIR(+AHD) epidemic model.
%  
%  Reference:
%  A. Cunha Jr , D. A. W. Barton, and T. G. Ritto
%  Uncertainty  quantification  in  epidemic  models  via
%  cross-entropy approximate Bayesian computation, 2022
% -----------------------------------------------------------------
%  programmers: Americo Cunha Jr (UERJ)
%               David A. W. Barton (Univ. Bristol)
%               Thiago G. Ritto (UFRJ)
%
%  last updated: May 26, 2022
% -----------------------------------------------------------------


clc
clear
close all

% program execution start time
% ----------------------------------------------------------------
timeStart = tic();
% ----------------------------------------------------------------


% program header
% ----------------------------------------------------------------
disp(' ---------------------------------------------------------- ')
disp(' Approximate Bayesian Computation (ABC)                     ')
disp(' for UQ in a SEIR(+AHD) model                               ')
disp('                                                            ')
disp(' by                                                         ')
disp(' Americo Cunha Jr (UERJ)                                    ')
disp(' David A. W. Barton (Univ. Briston)                         ')
disp(' Thiago G. Ritto (UFRJ)                                     ')
disp(' ---------------------------------------------------------- ')
% ----------------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'ABC_SEIRpAHD';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');

% random number generator (fix the seed for reproducibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);
% -----------------------------------------------------------


% load epidemic data
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- loading epidemic data --- ');
disp(' ');
disp('    ... ');
disp(' ');

load('COVID19_Data_RJ_2020_01_01_to_2021_01_01.mat')

% data range of dates
DateStart = datenum('01-01-2020');
DateEnd   = datenum('01-01-2021');

% training data range of dates
DateStart_train = datenum('05-01-2020');
DateEnd_train   = datenum('05-31-2020');

% validation data range of dates
DateStart_valid = DateEnd_train + 1;
DateEnd_valid   = DateStart_valid + 30;

% indices for training data beginning/end
Day0_train   = DateStart_train - DateStart;
DayEnd_train = DateEnd_train   - DateStart;

% indices for validation data beginning/end
Day0_valid   = DateStart_valid - DateStart;
DayEnd_valid = DateEnd_valid   - DateStart;

% time series of hospitalizations
DataH = Data_Hospitalized;

% time series of total deaths
DataD = cumsum(Data_NewDeaths);

% time series of hospitalizations (training)
DataH_train = DataH(Day0_train:DayEnd_train);

% time series of total deaths (training)
DataD_train = DataD(Day0_train:DayEnd_train);

% time series of hospitalizations (validation)
DataH_valid = DataH(Day0_valid:DayEnd_valid);

% time series of total deaths (validation)
DataD_valid = DataD(Day0_valid:DayEnd_valid);

% training data vector
data_train  = [DataH_train DataD_train];

% training data size / number of quantities of interest
[Ntrain,Nqoi] = size(data_train);

% validation data size 
Nvalid = length(DataH_valid);

% weights parameters
weigths = [0.75; 0.25];

toc
% -----------------------------------------------------------


% define model hyperparameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining model hyperparameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% initial time (days)
t0 = 1;

% final time (days)
t1 = t0 + DayEnd_valid - Day0_train;

% time step (days)
dt = 1;

% interval of analysis
tspan = t0:dt:t1;

% interval of analysis for a virgin population
tspan0 = 1:1:300;

% number of time steps
Ndt = length(tspan);

toc
% -----------------------------------------------------------


% define nominal model parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining model parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% transmission rate (days^-1)
beta     = 1/7;
beta_min = 1/14;
beta_max = 1/2;

% latent rate (days^-1)
alpha     = 1/5;
alpha_min = 1/10;
alpha_max = 1/2;

% symptomatic fraction (dimensionless)
% (0 <= fE <= 1)
fE     = 0.8;
fE_min = 0.7;
fE_max = 0.9;

% recovery rate (days^-1)
gamma     = 1/14;
gamma_min = 1/21;
gamma_max = 1/7;

% hospitalization rate (days^-1)
rho     = 0.01*(1/7);
rho_min = 0.01*(1/21);
rho_max = 0.01*(1/1);
% rho     = (1/6);
% rho_min = (1/21);
% rho_max = (1/5);

% death rate (days^-1)
delta     = 0.001*(1/14);
delta_min = 0.001*(1/21);
delta_max = 0.010*(1/1);
% delta     = (1/12);
% delta_min = (1/21);
% delta_max = (1/1);

% asymptomatic mortality-factor (dimensionless)
% (0 <= kappaA <= 1)
kappaA     = 0.0010;
kappaA_min = 0.0005;
kappaA_max = 0.0050;
% kappaA     = 0.75;
% kappaA_min = 0.5;
% kappaA_max = 1.0;

% hospitalization mortality-factor (dimensionless)
% (0 <= kappaH <= 1)
kappaH     = 0.05;
kappaH_min = 0.01;
kappaH_max = 0.10;
% kappaH     = 1.2;
% kappaH_min = 0.8;
% kappaH_max = 2.0;

% hospitalization infectivity-factor (dimensionless)
% (0 <= epsilonH <= 1)
epsilonH     = 0.2;
epsilonH_min = 0.1;
epsilonH_max = 0.5;

% parameters vector
param = [beta alpha fE gamma rho delta ...
         kappaA kappaH epsilonH];

% initial conditions for a virgin population
N0 = 5.5e6;                % initial population   (number of individuals)
D0 = 0;                    % initial deaths       (number of individuals)
H0 = 0;                    % initial hospitalized (number of individuals)
A0 = 0;                    % initial asymptomatic (number of individuals)
R0 = 0;                    % initial recovered    (number of individuals)
I0 = 0;                    % initial infectious   (number of individuals)
E0 = 1;                    % initial exposed      (number of individuals)
S0 = N0-D0-H0-A0-R0-I0-E0; % initial susceptible  (number of individuals)

% initial condition for a virgin population
IC0 = [S0 E0 I0 R0 A0 H0 D0 N0];

% initial condition compatible with initial deaths
IC = GetState_SEIRpAHD(param,tspan0,IC0,DataH_train(1),DataD_train(1));

% model function
fun = @(x) MyModel1(x,tspan,IC);

toc
% -----------------------------------------------------------



% ABC computation
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- ABC computation --- ');
disp(' ');
disp('    ... ');
disp(' ');

% parameters bounds
lb = [    beta_min; ...
         alpha_min; ...
            fE_min; ...
         gamma_min; ...
           rho_min; ...
         delta_min; ...
        kappaA_min; ...
        kappaH_min; ...
      epsilonH_min];
ub = [    beta_max; ...
         alpha_max; ...
            fE_max; ...
         gamma_max; ...
           rho_max; ...
         delta_max; ...
        kappaA_max; ...
        kappaH_max; ...
      epsilonH_max];

% parameters low-order statistics
mu    = [beta alpha fE gamma rho delta kappaA kappaH epsilonH];
sigma = (ub-lb)/sqrt(12);

% number of ABC samples
ABCobj.Ns= 2000;

% error tolerance for ABC
ABCobj.tol = 0.1;

% ABC weights parameters
ABCobj.weigths = weigths;

% prior for ABC
ABCobj.prior = 'TruncGaussian';
%ABCobj.prior = 'LogNormal';
%ABCobj.prior = 'Gamma';
%ABCobj.prior = 'Uniform';

% ABC algorithm
[x_best,y_best,ABCobj] = ABC(fun,data_train,lb,ub,mu,sigma,ABCobj);

toc
% -----------------------------------------------------------



% compute statistics
% -----------------------------------------------------------
tic
disp(' ')
disp(' --- computing statistics --- ');
disp(' ');
disp('    ... ');
disp(' ');

% quantities of interest best estimation
QoI_H_best = y_best(:,1);
QoI_D_best = y_best(:,2);

% quantities of interest accepted samples
QoI_H_samples = ABCobj.y_accept(1:Ndt    ,:)';
QoI_D_samples = ABCobj.y_accept(Ndt+1:end,:)';

% model parameters accepted samples
param_samples = ABCobj.x_accept;

% confidence probability (percentual)
Pc = 95;

% upper percentil
r_plus = 0.5*(100 + Pc);

% lower percentil
r_minus = 0.5*(100 - Pc);

% confidence bands lower bounds
QoI_H_low = prctile(QoI_H_samples,r_minus);
QoI_D_low = prctile(QoI_D_samples,r_minus);

% confidence bands upper bounds
QoI_H_upp = prctile(QoI_H_samples,r_plus);
QoI_D_upp = prctile(QoI_D_samples,r_plus);

% median
QoI_H_median = median(QoI_H_samples);
QoI_D_median = median(QoI_D_samples);

toc
% -----------------------------------------------------------


% integration of the system dynamics
%---------------------------------------------------------------
tic
disp(' '); 
disp(' --- integration of the system dynamics --- ');
disp(' ');
disp('    ... ');
disp(' ');

% ODE solver Runge-Kutta45    
[time,y] = ode45(@(t,y)rhs_SEIRpAHD(t,y,x_best),tspan,IC);

% define time series
S = y(:,1);  % susceptible         (number of individuals)
E = y(:,2);  % exposed             (number of individuals)
I = y(:,3);  % infected            (number of individuals)
R = y(:,4);  % recovered           (number of individuals)
A = y(:,5);  % asymptomatic        (number of individuals)
H = y(:,6);  % hospitalized        (number of individuals)
D = y(:,7);  % deaths              (number of individuals)
N = y(:,8);  % current population  (number of individuals)

toc
% -----------------------------------------------------------


% save simulation results
% -----------------------------------------------------------
tic
disp(' ')
disp(' --- saving simulation results --- ');
disp(' ');
disp('    ... ');
disp(' ');

save([num2str(case_name),'.mat']);

toc
% -----------------------------------------------------------



% post-processing
% -----------------------------------------------------------
tic
disp(' ')
disp(' --- post-processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


% ..........................................................
disp(' ');
disp(' ...........................');
disp('  ABC Identified Parameters ')
disp(' ...........................');
disp(['  beta     = ',num2str(x_best(1))])
disp(['  alpha    = ',num2str(x_best(2))])
disp(['  fE       = ',num2str(x_best(3))])
disp(['  gamma    = ',num2str(x_best(4))])
disp(['  rho      = ',num2str(x_best(5))])
disp(['  delta    = ',num2str(x_best(6))])
disp(['  kappaA   = ',num2str(x_best(7))])
disp(['  kappaH   = ',num2str(x_best(8))])
disp(['  epsilonH = ',num2str(x_best(9))])
disp(' ..........................');
disp(' ');
disp(' ..........................');
disp('  ABC algorithm statistics ')
disp(' ..........................');
disp(['  number of samples = ',num2str(ABCobj.Ns)])
disp(['  accepted  samples = ',num2str(ABCobj.accept_counter)])
disp(['  acceptance rate   = ',num2str(ABCobj.accept_rate)])
disp(['  error tolerance   = ',num2str(ABCobj.tol)])
disp(['  minimum error     = ',num2str(ABCobj.J_min)])
disp(' ..........................');
disp(' ');
% ..........................................................


% custom colors
% ..........................................................
yellow = [255 204  0]/256;
orange = [256 128  0]/256;
brown  = [101  33 33]/256;
%brownT = [255 235 205]/256;
brownT = [188 143 143]/256;
blackT = [192 192 192]/256;
% ..........................................................


% legend labels
% ..........................................................
graphobj.leg1 = 'Data: training';
graphobj.leg2 = 'Data: validation';
graphobj.leg3 = 'CE Optimal Fit';
graphobj.leg4 = 'ABC Best Fit';
graphobj.leg5 = 'ABC Median';
graphobj.leg6 = '95% envelope';
graphobj.leg7 = 'ABC samples';
% ..........................................................


% figs title
% ..........................................................
fig_title = 'SEIR(+AHD) model ';
% ..........................................................


% adjust time vector for date format
% ..........................................................
time       = linspace(DateStart_train,DateEnd_valid,Ndt);
time_train = linspace(DateStart_train,DateEnd_train,Ntrain);
time_valid = linspace(DateStart_valid,DateEnd_valid,Nvalid);
% ..........................................................


% ..........................................................
% plot the quantities of interest (hospitalizations)
% ..........................................................
graphobj.gname = [num2str(case_name),'__UQ_H'];
graphobj.gtitle = '';
graphobj.xmin   = DateStart_train;
graphobj.xmax   = DateEnd_valid;
graphobj.ymin   = 0;
graphobj.ymax   = 1e3;
graphobj.xlab   = [];
graphobj.ylab   = 'hospitalized individuals';
graphobj.flag   = 'eps';
graphobj.color  = brown;
graphobj.colorT = brownT;

fig_UQ_H = graph_QoI_UQ(time_train,DataH_train,...
                        time_valid,DataH_valid,...
                                 QoI_H_samples,...
                              time,QoI_H_opt,...
                                   QoI_H_best,...
                                   QoI_H_median,...
                                   QoI_H_low,...
                                   QoI_H_upp,...
                                   graphobj);
% ..........................................................


% ..........................................................
% plot the quantities of interest (total deaths)
% ..........................................................
graphobj.gname = [num2str(case_name),'__UQ_D'];
graphobj.gtitle = '';
graphobj.xmin   = DateStart_train;
graphobj.xmax   = DateEnd_valid;
graphobj.ymin   = 0;
graphobj.ymax   = 30e3;
graphobj.xlab   = [];
graphobj.ylab   = 'total deaths';
graphobj.flag   = 'eps';
graphobj.color  = 'k';
graphobj.colorT = blackT;

fig_UQ_D = graph_QoI_UQ(time_train,DataD_train,...
                        time_valid,DataD_valid,...
                                 QoI_D_samples,...
                              time,QoI_D_opt,...
                                   QoI_D_best,...
                                   QoI_D_median,...
                                   QoI_D_low,...
                                   QoI_D_upp,...
                                   graphobj);
% ..........................................................


% ..........................................................
% plot parameters samples from ABC
% ..........................................................
graphobj.gname = [num2str(case_name),'__Samples'];
graphobj.gtitle = '';
graphobj.xmin   = 'auto';
graphobj.xmax   = 'auto';
graphobj.ymin   = 'auto';
graphobj.ymax   = 'auto';
graphobj.xlab   = [];
graphobj.ylab   = [];
graphobj.flag   = 'eps';

%fig_ABC_samples = graph_SampleMatrix(ABCobj,graphobj);
% ..........................................................

% ..........................................................
% plot all populations of SEIR(+AHD) model (log scale)
% ..........................................................
graphobj.gname  = [num2str(case_name),'__time_series_log'];
graphobj.labelS = 'Suceptibles';
graphobj.labelE = 'Exposed';
graphobj.labelI = 'Infected';
graphobj.labelR = 'Recovered';
graphobj.labelA = 'Asymptomatic';
graphobj.labelH = 'Hospitalized';
graphobj.labelD = 'Deaths';
graphobj.labelN = 'Population';
graphobj.gtitle = '';
graphobj.xmin   = DateStart_train;
graphobj.xmax   = DateEnd_valid;
graphobj.ymin   = 1;
graphobj.ymax   = 'auto';
graphobj.xlab   = [];
graphobj.ylab   = 'number of individuals';
graphobj.flag   = 'eps';

fig_SEIRpAHDlog_opt = graph_SEIRpAHDlog_opt(time,S,E,I,R,A,H,D,N,graphobj);
% ..........................................................

toc
% -----------------------------------------------------------


% program execution time
% -----------------------------------------------------------
disp(' ');
disp(' -----------------------------');
disp('            THE END!          ');
disp(' -----------------------------');
disp('  Total execution time:       ');
disp(['  ',num2str(toc(timeStart)),' seconds']);
disp(' -----------------------------');
% -----------------------------------------------------------
