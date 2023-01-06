% -----------------------------------------------------------------
% GetState_SEIRpAHD.m
% -----------------------------------------------------------------
%  This function computes a dynamic state that is compatible with
%  a given death reference for a SEIR(+AHD) epidemic model.
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
%    beta     = transmission rate                  (days^-1)
%    alpha    = latent rate                        (days^-1)
%    fE       = symptomatic fraction               (dimensionless)
%    gamma    = recovery rate                      (days^-1)
%    rho      = hospitalization rate               (days^-1)
%    delta    = death rate                         (days^-1)
%    kappaH   = hospitalization mortality-factor   (dimensionless)
%    epsilonH = hospitalization infectivity-factor (dimensionless)
%  
%  This code was adapted from:
%  EPIDEMIC - Epidemiology Educational Code
%  www.EpidemicCode.org
%  
%  Reference:
%  A. Cunha Jr , T. G. Ritto, and D. A. W. Barton,
%  Uncertainty  quantification  in  epidemic  models  via
%  cross-entropy approximate Bayesian computation
%  Preprint, 2022
% -----------------------------------------------------------------
%  programmers: Americo Cunha Jr (UERJ)
%               David A. W. Barton (Univ. Bristol)
%               Thiago G. Ritto (UFRJ)
%
%  last updated: Dec 20, 2022
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [IC,IC_H,idx_H0,IC_D,idx_D0] = ...
                  GetState_SEIRpAHD(param,tspan,IC0,H_ref,D_ref)

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
N = y(:,8);  % population          (number of individuals)

% under notification factor
U = 1.0;

% dynamic state near the reference H_ref
[H0dif,idx_H0] = min(abs(H-H_ref));

N0  = IC0(8);
%D0  = D(idx_H0,:)*U;
D0  = D_ref;
H0  = H_ref;
A0  = A(idx_H0,:)*U;
R0  = R(idx_H0,:)*U;
I0  = I(idx_H0,:)*U;
E0  = E(idx_H0,:)*U;
S0  = N0-D0-H0-A0-R0-I0-E0;
IC_H = [S0 E0 I0 R0 A0 H0 D0 N0];

% dynamic state near the reference D_ref
[D0dif,idx_D0] = min(abs(D-D_ref));

N0  = IC0(8);
D0  = D_ref;
H0  = H_ref;
%H0  = H(idx_D0,:)*U;
A0  = A(idx_D0,:)*U;
R0  = R(idx_D0,:)*U;
I0  = I(idx_D0,:)*U;
E0  = E(idx_D0,:)*U;
S0  = N0-D0-H0-A0-R0-I0-E0;
IC_D = [S0 E0 I0 R0 A0 H0 D0 N0];

% weigthed dynamic state
IC = 0.5*(IC_H + IC_D);

end
% -----------------------------------------------------------------

