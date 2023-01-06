% -----------------------------------------------------------------
%  MyModel2.m
% -----------------------------------------------------------------
%  This function defines the model function. It must be customized 
%  according to your particular application.
%  
%  In this case it is adjusted to a SEIR(+AHD)beta epidemic model.
%
%  Input:
%  x     - model parameters vector
%  tspan - temporal mesh vector
%  IC    - initial conditions vector
%  
%  Output:
%  F - matrix with the quantities of interest time series
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
function F = MyModel2(x,tspan,IC)

    % preallocate memory for the model response
    F = zeros(length(tspan),2);

	% ODE solver Runge-Kutta45
	[time,ymodel] = ode45(@(t,y)rhs_SEIRpAHDbeta(t,y,x),tspan,IC);
    
    % model response
    F(:,1) = ymodel(:,6);
    F(:,2) = ymodel(:,7);
end
% -----------------------------------------------------------------