% -----------------------------------------------------------------
%  Main_COVID19RJ_Data_plot.m
% -----------------------------------------------------------------
%  This program plots COVID-19 data for Rio de Janeiro city.
%  
%  Reference:
%  A. Cunha Jr , D. A. W. Barton, and T. G. Ritto
%  Uncertainty  quantification  in  epidemic  models  via
%  cross-entropy approximate Bayesian computation, 2023
% -----------------------------------------------------------------
%  programmers: Americo Cunha Jr (UERJ)
%               David A. W. Barton (Univ. Bristol)
%               Thiago G. Ritto (UFRJ)
%
%  last updated: Dec 20, 2022
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
disp(' COVID-19 data for Rio de Janeiro city                      ')
disp('                                                            ')
disp(' by                                                         ')
disp(' Americo Cunha Jr (UERJ)                                    ')
disp(' David A. W. Barton (Univ. Briston)                         ')
disp(' Thiago G. Ritto (UFRJ)                                     ')
disp(' ---------------------------------------------------------- ')
% ----------------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'COVID19RJ_Data';

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

% data vectors
DataH = Data_Hospitalized;
DataD = cumsum(Data_NewDeaths);
data  = [DataH DataD];

% data time series size / number of quantities of interest
[Ndata,Nqoi] = size(data);

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
t1 = t0 + DateEnd - DateStart;

% time step (days)
dt = 1;

% interval of analysis
tspan = t0:dt:t1;

% temporal mesh
time = tspan';

% number of time steps
Ndt = length(tspan);

% number of time steps for a single day/week
Nday  = round(1/dt);

% cumulative number of hospitalizations
DataHcum = cumsum(DataH);

% news deaths per day (number of individuals/day)
DataDnew          = zeros(Ndt,1);
DataDnew(1:end-1) = DataD(2:Nday:end)-DataD(1:Nday:Ndt-Nday);
DataDnew(end)     = DataDnew(end-1);

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

% custom colors
% ..........................................................
brown  = [101  33 33]/256;
% ..........................................................


% adjust time vector for date format
% ..........................................................
time = linspace(DateStart,DateEnd,length(time));
% ..........................................................


% ..........................................................
% plot H data
% ..........................................................
graphobj.gname = [num2str(case_name),'_cum'];
graphobj.gtitle = '';
graphobj.xmin   = DateStart;
graphobj.xmax   = DateEnd;
graphobj.ymin   = 0;
graphobj.ymax   = 'auto';
graphobj.xlab   = [];
graphobj.ylab   = 'number of individuals';
graphobj.flag   = 'eps';
graphobj.color1 = brown;
graphobj.color2 = 'k';
graphobj.leg1   = 'total deaths';

fig_Data_cum = graph_Data1(time,DataD,graphobj);
% ..........................................................


% ..........................................................
% plot D data
% ..........................................................
graphobj.gname = [num2str(case_name),'_new'];
graphobj.gtitle = '';
graphobj.xmin   = DateStart;
graphobj.xmax   = DateEnd;
graphobj.ymin   = 0;
graphobj.ymax   = 1e3;
graphobj.xlab   = [];
graphobj.ylab   = 'number of individuals';
graphobj.flag   = 'eps';
graphobj.color1 = brown;
graphobj.color2 = 'k';
graphobj.leg1   = 'hospitalized';
graphobj.leg2   = 'new deaths';

fig_Data_new = graph_Data2(time,DataH,DataDnew,graphobj);
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