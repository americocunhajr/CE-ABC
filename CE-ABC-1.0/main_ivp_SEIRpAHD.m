% -----------------------------------------------------------------
%  main_ivp_SEIRpAHD.m
% -----------------------------------------------------------------
%  This program simulates the nonlinear dynamics of an 
%  SEIR(+AHD) epidemic model.
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
%  last updated: April 22, 2022
% -----------------------------------------------------------------


clc
clear
close all


% program header
% ----------------------------------------------------------------
disp(' -------------------------------------------------------- ')
disp(' SEIR(+AHD) epidemic model                                ')
disp(' (dynamic integration)                                    ')
disp('                                                          ')
disp(' by                                                       ')
disp(' Thiago G. Ritto (UFRJ)                                   ')
disp(' Americo Cunha Jr (UERJ)                                  ')
disp(' David A. W. Barton (Univ. Briston)                       ')
disp(' -------------------------------------------------------- ')
% ----------------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'SEIRpAHD_ivp';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
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
DateTrainStart = datenum('05-01-2020');
DateTrainEnd   = datenum('01-01-2021');

% indices for training data beginning/end
Day0   = DateTrainStart - DateStart;
DayEnd = DateTrainEnd   - DateStart;

% data vectors
DataH = Data_Hospitalized(Day0:DayEnd);
DataD = cumsum(Data_NewDeaths(Day0:DayEnd));
data  = [DataH DataD];

% data time series size / number of quantities of interest
[Ndata,Nqoi] = size(data);

toc
% -----------------------------------------------------------



% define model parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining model parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% transmission rate (days^-1)
beta = 1/7;

% latent rate (days^-1)
alpha = 1/5;

% symptomatic fraction (dimensionless)
% (0 <= fE <= 1)
fE = 0.8;

% recovery rate (days^-1)
gamma = 1/14;

% hospitalization rate (days^-1)
%rho = (1/100)*(1/15);
rho = (1/100)*(1/6);

% death rate (days^-1)
%delta = (1/1000)*(1/7);
delta = (1/1000)*(1/15);

% asymptomatic mortality-factor (dimensionless)
% (0 <= kappaA <= 1)
kappaA = 0.001;

% hospitalization mortality-factor (dimensionless)
% (0 <= kappaH <= 1)
kappaH = 0.05;

% hospitalization infectivity-factor (dimensionless)
% (0 <= epsilonH <= 1)
epsilonH = 0.2;

% parameters vector
param = [beta alpha fE gamma rho delta kappaA kappaH epsilonH];

% initial conditions
N0 = 5.5e6;                % initial population   (number of individuals)
D0 = 0;                    % initial deaths       (number of individuals)
H0 = 0;                    % initial hospitalized (number of individuals)
A0 = 0;                    % initial asymptomatic (number of individuals)
R0 = 0;                    % initial recovered    (number of individuals)
I0 = 0;                    % initial infectious   (number of individuals)
E0 = 1;                    % initial exposed      (number of individuals)
S0 = N0-D0-H0-A0-R0-I0-E0; % initial susceptible  (number of individuals)

% initial conditions vector
IC0 = [S0 E0 I0 R0 A0 H0 D0 N0];

% basic reproduction number
R_nought = beta/gamma;

% mean duration in compartment I
DI = gamma + rho + delta;

% mean duration in compartment A
DA = gamma + delta;

% mean duration in compartment H
DH = gamma + kappaH*delta;

% control reproduction number
R_control = fE*beta/DI + (1-fE)*beta/DA;

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

% final time   (days)
t1 = 2*365;

% time step (days)
dt = 1;

% interval of analysis
tspan = t0:dt:t1;

% number of time steps
Ndt = length(tspan);

% number of time steps for a single day/week
Nday = round(1/dt);

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
[time, y] = ode45(@(t,y)rhs_SEIRpAHD(t,y,param),tspan,IC0);

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


% computing the quantities of interest 
%---------------------------------------------------------------
tic
disp(' '); 
disp(' --- computing the quantities of interest --- ');
disp(' ');
disp('    ... ');
disp(' ');

% cumulative number of hospitalizations
Hcum = Map_New2Cum([H(1); rho*I(1:end-1)]);

% news deaths per day (number of individuals/day)
Dnew = Map_Cum2New(D);

% states of reference
[IC,IC_H,idx_H0,IC_D,idx_D0] = ...
    state_SEIRpAHD(param,1:1:365,IC0,DataH(1),DataD(1));

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
% plot all populations of SEIR(+AHD) model
% ..........................................................
graphobj.gname  = [num2str(case_name),'__time_series'];
graphobj.labelS = 'Suceptibles';
graphobj.labelE = 'Exposed';
graphobj.labelI = 'Infected';
graphobj.labelR = 'Recovered';
graphobj.labelA = 'Asymptomatic';
graphobj.labelH = 'Hospitalized';
graphobj.labelD = 'Deaths';
graphobj.labelN = 'Population';
graphobj.gtitle = '';
graphobj.xmin   = t0;
graphobj.xmax   = t1;
graphobj.ymin   = 0;
graphobj.ymax   = 6e6;
graphobj.xlab   = 'time (days)';
graphobj.ylab   = 'number of individuals';
graphobj.flag   = 'eps';

fig_SEIRpAHD = graph_SEIRpAHD(time,S,E,I,R,A,H,D,N,graphobj);
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
graphobj.xmin   = t0;
graphobj.xmax   = t1;
graphobj.ymin   = 1;
graphobj.ymax   = 10e6;
graphobj.xlab   = 'time (days)';
graphobj.ylab   = 'number of individuals';
graphobj.flag   = 'eps';

fig_SEIRpAHDlog = graph_SEIRpAHDlog(time,S,E,I,R,A,H,D,N,graphobj);
% ..........................................................


% ..........................................................
% plot the quantities of interest (cumulative)
% ..........................................................
graphobj.gname     = [num2str(case_name),'__QoIcum'];
graphobj.labelCumH = 'Cumulative Hospitalizations';
graphobj.labelD    = 'Deaths';
graphobj.gtitle    = '';
graphobj.xmin      = t0;
graphobj.xmax      = t1;
graphobj.ymin      = 0;
graphobj.ymax      = 50e3;
graphobj.xlab      = 'time (days)';
graphobj.ylab      = 'number of individuals';
graphobj.flag      = 'eps';

fig_QoIcum = graph_QoIcum(time,Hcum,D,graphobj);
% ..........................................................


% ..........................................................
% plot the quantities of interest (new number)
% ..........................................................
graphobj.gname     = [num2str(case_name),'__QoInew'];
graphobj.labelH    = 'Hospitalized';
graphobj.labelNewD = 'New Deaths';
graphobj.gtitle    = '';
graphobj.xmin      = t0;
graphobj.xmax      = t1;
graphobj.ymin      = 0;
graphobj.ymax      = 7e3;
graphobj.xlab      = 'time (days)';
graphobj.ylab      = 'number of individuals';
graphobj.flag      = 'eps';

fig_QoInew = graph_QoInew(time,H,Dnew,graphobj);
% ..........................................................


% ..........................................................
% plot all populations of SEIR(+AHD) model with ref states
% ..........................................................
graphobj.gname      = [num2str(case_name),'__ref_states'];
graphobj.labelDataH = 'Ref. State: Hospitalized';
graphobj.labelDataD = 'Ref. State: Deaths';
graphobj.gtitle     = '';
graphobj.xmin       = t0;
graphobj.xmax       = t1;
graphobj.ymin       = 1;
graphobj.ymax       = 10e6;
graphobj.xlab       = 'time (days)';
graphobj.ylab       = 'number of individuals';
graphobj.flag       = 'eps';

fig_SEIRpAHDstate = ...
 graph_SEIRpAHDstate(time,t0,t1,S,E,I,R,A,H,D,N,...
                     DataH(1),DataD(1),idx_H0,idx_D0,graphobj);
% ..........................................................

toc
% -----------------------------------------------------------

