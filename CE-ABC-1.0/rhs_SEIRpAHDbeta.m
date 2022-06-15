% -----------------------------------------------------------------
% rhs_SEIRpAHDbeta.m
% -----------------------------------------------------------------
%  This function defines the system of ODEs for the SEIR(+AHD)beta
%  epidemic model.
%
%  The dynamic state coordinates are:
%    S = susceptibles       (number of individuals)
%    E = exposed            (number of individuals)
%    I = infected           (number of individuals)
%    R = recovered          (number of individuals)
%    A = asymptomatic       (number of individuals)
%    H = hospitalized       (number of individuals)
%    D = deaths             (number of individuals)
%    N = living population  (number of individuals)
%
%  The epidemic model parameters are:
%    beta0    = initial transmission rate          (days^-1)
%    alpha    = latent rate                        (days^-1)
%    fE       = symptomatic fraction               (dimensionless)
%    gamma    = recovery rate                      (days^-1)
%    rho      = hospitalization rate               (days^-1)
%    delta    = death rate                         (days^-1)
%    kappaA   = asymptomatic mortality-factor      (dimensionless)
%    kappaH   = hospitalization mortality-factor   (dimensionless)
%    epsilonH = hospitalization infectivity-factor (dimensionless)
%    beta_inf = final transmission rate            (days^-1)
%    eta      = rate of change for beta            (days^-1)
%    tau_beta = tranmission rate inflexion time    (day    )
%  
%  This code was adapted from:
%  EPIDEMIC - Epidemiology Educational Code
%  www.EpidemicCode.org
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
%  last update: March 17, 2022
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function dydt = rhs_SEIRpAHDbeta(t,y,param)

% number of parameters
Nparam = length(param);

% number of jumpes between waves
Njumps = round((Nparam-9)/3);

% SEIR(+AHD) dynamic model parameters:
%   beta0    - initial transmission rate          (days^-1)
%   alpha    - latent rate                        (days^-1)
%   fE       - symptomatic fraction               (dimensionless)
%   gamma    - recovery rate                      (days^-1)
%   rho      - hospitalization rate               (days^-1)
%   delta    - mortality rate                     (days^-1)
%   kappaA   - asymptomatic mortality-factor      (dimensionless)
%   kappaH   - hospitalization mortality-factor   (dimensionless)
%   epsilonH - hospitalization infectivity-factor (dimensionless)
%   beta_inf - final transmission rate            (days^-1)
%   eta      - rate of change for beta            (days^-1)
%   tau_beta - tranmission rate inflexion time    (day    )

beta0    = param(1);
alpha    = param(2);
fE       = param(3);
gamma    = param(4);
rho      = param(5);
delta    = param(6);
kappaA   = param(7);
kappaH   = param(8);
epsilonH = param(9);
beta_inf = param(10:3:Nparam);
eta      = param(11:3:Nparam);
tau_beta = param(12:3:Nparam);

beta = beta0;
for i=1:Njumps
    beta  = beta + 0.5*(beta_inf(i)-beta0)*(1+tanh(0.5*eta(i)*(t-tau_beta(i))));
    beta0 = beta_inf(i);
end

% SEIR(+AHD) dynamic model equations:
% 
%    y = [S E I R A H D N]                         is the state vector
% dydt = [dSdt dEdt dIdt dRdt dAdt dHdt dDdt dNdt] is the evolution law
% 
% dSdt - rate of susceptible          (number of individuals/days)
% dEdt - rate of exposed              (number of individuals/days)
% dIdt - rate of infected             (number of individuals/days)
% dRdt - rate of recovered            (number of individuals/days)
% dAdt - rate of asymptomatic         (number of individuals/days)
% dHdt - rate of hospitalized         (number of individuals/days)
% dDdt - rate of deaths               (number of individuals/days)
% dNdt - rate of population variation (number of individuals/days)

[S E I R A H D N] = deal(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8));

% model ODEs
dSdt = - beta*S.*(I+A+epsilonH*H)./N;
dEdt =   beta*S.*(I+A+epsilonH*H)./N - alpha*E;
dIdt = fE*alpha*E - (rho+delta+gamma)*I;
dRdt = gamma*(I+A+H);
dAdt = (1-fE)*alpha*E - (kappaA*delta+gamma)*A;
dHdt = rho*I - (gamma+kappaH*delta)*H;
dDdt = delta*(I+kappaA*A+kappaH*H);
dNdt = -dDdt;

% system of ODEs
dydt = [dSdt; dEdt; dIdt; dRdt; dAdt; dHdt; dDdt; dNdt];

end
% -----------------------------------------------------------------