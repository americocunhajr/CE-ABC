% -----------------------------------------------------------------
%  MyCostFunc.m
% -----------------------------------------------------------------
%  This function defines the cost function to be called by
%  CE algorithm. It must be customized according to your
%  particular application.
%  
%  In this case it is adjusted to a compute the discrepancy
%  between SEIR(+AHD) epidemic model prediction and a given
%  set of epidemic data.
%
%  Input:
%  x       - design variables vector
%  fun     - model prediction function
%  data    - data vector
%  weigths - weigths vector
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
%  last update: April 25, 2022
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function J = MyCostFunc(x,fun,data,weigths)

    % input dimensions
    [Nvars,Ns] = size(x);
    
    % number of data points and outputs
    [Ndata,Nqoi] = size(data);

    % preallocate memory for the cost function
    J = zeros(1,Ns);
    
    % squared norm of data vector
    norm_data_pow2 = data'*data;
    
    % loop over samples to compute the cost function
    for n = 1:Ns
        
        % model response
        y = fun(x(:,n));
        
        % difference between model and data
        delta_y = data - y(1:Ndata,1:Nqoi);

        % CE cost function
        J(n) = sum(weigths.*diag(((delta_y')*delta_y)./norm_data_pow2));
    end
end
% -----------------------------------------------------------------