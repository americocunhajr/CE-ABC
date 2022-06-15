% -----------------------------------------------------------------
%  CEopt.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last updated: April 23, 2022
% -----------------------------------------------------------------
%  This routine employs the Cross-entropy method to solve the 
%  following optimization problem
%  
%         X_opt = arg min F(x)
%             lb <= x <= up
%  
%  where F: R^{Nvars} --> R is a given scalar function.
%  
%  The goal here is to minimize a scalar objective function defined
%  in a known rectangular domain (feasible region). For constrained
%  optimization problems, formulations based on penalty method (or 
%  augmented Lagrangian method) is expected, the user must provide
%  an unconstrained objective function, which incorporates in its
%  definition of all the equality constraints of the original problem.
%  
%  The employed algorithm samples the feasible region using a
%  truncated Gaussian distribution and update its parameters, with
%  aid of an elite set (defined by the better samples), seeking to
%  transform (in the limit) this Gaussian into a Dirac distribution
%  centered in the global optimum.
%  
%  Input:
%  fun    - objective function
%  x0     - initial mean
%  sigma0 - intial standard deviation
%  lb     - lower bound
%  ub     - upper bound
%  CEobj  - Cross-Entropy object struct
%  
%  Output:
%  X_opt - optimal point
%  F_opt - optimal value
%  CEobj - Cross-Entropy object struct
%  
%  Reference:
%  Reuven Y. Rubinstein, Dirk P. Kroese
%  The Cross-Entropy Method: A Unified Approach to Combinatorial 
%  Optimization, Monte-Carlo Simulation, and Machine Learning
%  Springer-Verlag, 2004.
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [X_opt,F_opt,CEobj] = CEopt(fun,x0,sigma0,lb,ub,CEobj)

    % check number of arguments
    if nargin < 5
        error('Too few inputs.')
    elseif nargin > 6
        error('Too many inputs.')
    end
    
   % check for consistency in lb and ub
    if size(lb) ~= size(ub)
        error('lb and ub must have the same dimensions')
    end
    if isempty(lb) || isempty(ub)
        error('lb and ub must be non-empty')
    end
    if any(isnan(lb)) || any(~isfinite(lb))
        error('lb must have finite values')
    end
    if any(isnan(ub)) || any(~isfinite(ub))
        error('ub must have finite values')
    end
    if any(lb > ub)
        error('lb < ub for all components')
    end
    
    % convert lb to a column vector (if necessary)
    [s1,s2] = size(lb);
    if isrow(lb)
        lb = lb(:);
    elseif s2 ~= 1
        error('lb must be a column vector')
    end
    
    % convert ub to a column vector (if necessary)
    [s1,s2] = size(ub);
    if isrow(ub)
        ub = ub(:);
    elseif s2 ~= 1 
        error('ub must be a column vector')
    end
    
    % define x0 (if necessary)
    if isempty(x0)
        x0 = (ub+lb)/2;
    end
    
    % define sigma0 (if necessary)
    if isempty(sigma0)
        sigma0 = (ub-lb)/sqrt(12);
    end
    
    % check for consistency in x0 and sigma0    
    if size(x0) ~= size(sigma0)
        error('x0 and sigma0 must have the same dimensions')
    end
    if any(isnan(x0)) || any(~isfinite(x0))
        error('x0 must have finite values')
    end
    if any(isnan(sigma0)) || any(~isfinite(sigma0))
        error('sigma0 must have finite values')
    end
    if any(x0 < lb) || any(x0 > ub)
        error('x0 must be in [lb,ub] interval')
    end
    if any(sigma0 <= 0.0)
        error('All components of sigma0 must be positive')
    end
    
    % convert x0 to a column vector (if necessary)
    [s1,s2] = size(x0);
    if isrow(x0)
        x0 = x0(:);
    elseif s2 ~= 1 
        error('x0 must be a column vector')
    end
    
    % convert sigma0 to a column vector (if necessary)
    [s1,s2] = size(sigma0);
    if isrow(sigma0)
        sigma0 = sigma0(:);
    elseif s2 ~= 1 
        error('sigma0 must be a column vector')
    end
    
    % number of variables
    Nvars = length(lb);
    
    % check for dimension consistency in x0
    if length(x0) ~= Nvars
        error('x0 must be a Nvars x 1 vector')
    end
    
    % check for dimension consistency in sigma0
    if length(sigma0) ~= Nvars
        error('sigma0 must be a Nvars x 1 vector')
    end
   
    % check if CEobj is declared
    if nargin == 5
        CEobj = [];
    end
        
    % maximum number of iterations
    if ~isfield(CEobj, 'maxiter')
        CEobj.maxiter = 100*Nvars;
    elseif CEobj.maxiter < 1
        error('maxiter must be a positive integer')
    end
    
    maxiter = CEobj.maxiter;

    % absolute tolerance
    if ~isfield(CEobj, 'atol')
        CEobj.atol = 1.0e-3;
    elseif CEobj.atol <= 0.0
        error('atol must be positive')
    end
    
    atol = CEobj.atol;

    % relative tolerance
    if ~isfield(CEobj, 'rtol')
        CEobj.rtol = 1.0e-3;
    elseif CEobj.rtol < 0.0
        error('rtol must be non-negative')
    end
    
    rtol = CEobj.rtol;

    % number of samples
    if ~isfield(CEobj, 'N')
        CEobj.N = 100;
    elseif CEobj.N < 1
        error('N must be greather than or equal to 1')
    end
    
    N = CEobj.N;

    % elite samples percentage
    if ~isfield(CEobj, 'rho')
        CEobj.rho = 0.05;
    elseif CEobj.rho <= 0 || CEobj.rho >= 1
        error('rho must be such that 0 < rho < 1')
    end
    
    rho = CEobj.rho;

    % smoothing parameter (0 < alpha <= 1)
    if ~isfield(CEobj, 'alpha')
        CEobj.alpha = 0.4;
    elseif CEobj.alpha < 0 || CEobj.alpha > 1
        error('alpha must be such that 0 <= alpha <= 1')
    end
    
    alpha = CEobj.alpha;

    % dynamic smoothing parameter
    if ~isfield(CEobj, 'beta')
        CEobj.beta = 0.4;
    elseif CEobj.beta <= 0
        error('beta must be non-negative')
    end
    
    beta = CEobj.beta;

    % dynamic smoothing parameter
    if ~isfield(CEobj, 'q')
        CEobj.q = 10;
    elseif CEobj.q <= 0
        error('q must be non-negative')
    end
    
    q = CEobj.q;

    % print on screen flag
    if ~isfield(CEobj, 'PRINT_ON')
        CEobj.PRINT_ON = 0;
    end
    
    PRINT_ON = CEobj.PRINT_ON;
        
    % number of samples in elite set
    N_elite = ceil(rho*N);

    % 1 - alpha variable
    one_minus_alpha = 1 - alpha;
    
    % initialize level counter
    t = 0;

    % preallocate memory for the ObjFunc value matrix
    CEobj.F = zeros(maxiter,1);
    
    % preallocate memory for the design variables matrix
    CEobj.x = zeros(maxiter,Nvars);

    % preallocate memory for the standard deviation matrix
    CEobj.sigma = zeros(maxiter,Nvars);

    % preallocate memory for the error weigthed vector
    CEobj.err_w = zeros(maxiter,Nvars);
    
    % preallocate memory for the error weighted norm
    CEobj.err_wrms = zeros(maxiter,1);
    
    % preallocate memory for design variables samples
    X = zeros(Nvars,N);
    
    % initialize ObjFunc minimum value
    F_opt = Inf;

    % initialize mean value
    mu = x0;

    % initialize standard deviation
    sigma = sigma0;
    
    % convergence indicator
    err_small = 0;
    
    % loop to sample the domain and update the distribution
    while err_small ~= 1 && t < maxiter

        % update level counter
        t = t + 1;

        % limit vectors for standard truncated Gaussian
        supp_l = ((lb - mu)./sigma)*ones(1,N);
        supp_u = ((ub - mu)./sigma)*ones(1,N);

        % generate samples from truncated Gaussian distribution
        for n=1:Nvars
            X(n,:) = mu(n) + trandn(supp_l(n,:),supp_u(n,:))*sigma(n);
        end

        % evaluate ObjFunc at the samples
        F = fun(X);

        % sort ObjFunc evaluations
        [F_sort,I_sort] = sort(F);
                
        % define elite sample set
        X_elite = X(:,I_sort(1:N_elite));
        
        % update mean value
        mu = mean(X_elite,2);
        
        % update standard deviation
        sigma = std(X_elite,0,2);

        % estimator for the ObjFunc minimum value
        F_min = F_sort(N_elite);

        % smoothing
        mu = alpha*mu + one_minus_alpha*x0;

        % update dynamics smoothing parameter
        beta_t = beta*(1 - (1-1/t)^q);

        % dynamic smoothing
        sigma = beta_t*sigma + (1-beta_t)*sigma0;
        
        % convergence metric
        [err_small,err_wrms] = ConvMetric(sigma,sigma0,atol,rtol);

        % update old mean
        x0 = mu;

        % update old standard deviation
        sigma0 = sigma;
        
        % update optimum
        if F_min < F_opt
           F_opt = F_min;
           X_opt = mu;
        end
        
        % save the algorithm records
        CEobj.F(t,1)        = F_min;
        CEobj.x(t,:)        = mu;
        CEobj.sigma(t,:)    = sigma;
        CEobj.err_wrms(t,:) = err_wrms;
        
        % print values on screen
        if PRINT_ON
            PrintScreen(t,F_min,err_wrms,mu);
        end
    end
    
    % print linebreak in the end
	if PRINT_ON
        fprintf('\n \n');
	end
    
	% delete empty entries from sampling records
    if t < maxiter
        CEobj.F        = CEobj.F(1:t,1);
        CEobj.x        = CEobj.x(1:t,:);
        CEobj.sigma    = CEobj.sigma(1:t,:);
        CEobj.err_wrms = CEobj.err_wrms(1:t,1);
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% trandn - truncated normal generator
% -----------------------------------------------------------------
% * efficient generator of a vector of length(l)=length(u)
% from the standard multivariate normal distribution,
% truncated over the region [l,u];
% infinite values for 'u' and 'l' are accepted;
% * Remark:
% If you wish to simulate a random variable
% 'Z' from the non-standard Gaussian N(m,s^2)
% conditional on l<Z<u, then first simulate
% X=trandn((l-m)/s,(u-m)/s) and set Z=m+s*X;
% 
% Reference: 
% Botev, Z. I. (2016). "The normal law under linear restrictions: 
% simulation and estimation via minimax tilting". Journal of the 
% Royal Statistical Society: Series B (Statistical Methodology). 
% doi:10.1111/rssb.12162
% -----------------------------------------------------------------
function x=trandn(l,u)
    l=l(:);u=u(:); % make 'l' and 'u' column vectors
    if length(l)~=length(u)
        error('Truncation limits have to be vectors of the same length')
    end
    x=nan(size(l));
    a=.66; % treshold for switching between methods
    % threshold can be tuned for maximum speed for each Matlab version
    % three cases to consider:
    % case 1: a<l<u
    I=l>a;
    if any(I)
        tl=l(I); tu=u(I); x(I)=ntail(tl,tu);
    end
    % case 2: l<u<-a
    J=u<-a;
    if any(J)
        tl=-u(J); tu=-l(J); x(J)=-ntail(tl,tu);
    end
    % case 3: otherwise use inverse transform or accept-reject
    I=~(I|J);
    if  any(I)
        tl=l(I); tu=u(I); x(I)=tn(tl,tu);
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% ntail - samples a column vector of length=length(l)=length(u)
% from the standard multivariate normal distribution,
% truncated over the region [l,u], where l>0 and
% l and u are column vectors;
% uses acceptance-rejection from Rayleigh distr. 
% similar to Marsaglia (1964);
% -----------------------------------------------------------------
function x=ntail(l,u)
    c=l.^2/2; n=length(l); f=expm1(c-u.^2/2);
    x=c-reallog(1+rand(n,1).*f); % sample using Rayleigh
    % keep list of rejected
    I=find(rand(n,1).^2.*x>c); d=length(I);
    while d>0 % while there are rejections
        cy=c(I); % find the thresholds of rejected
        y=cy-reallog(1+rand(d,1).*f(I));
        idx=rand(d,1).^2.*y<cy; % accepted
        x(I(idx))=y(idx); % store the accepted
        I=I(~idx); % remove accepted from list
        d=length(I); % number of rejected
    end
    x=sqrt(2*x); % this Rayleigh transform can be delayed till the end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% tn - samples a column vector of length=length(l)=length(u)
% from the standard multivariate normal distribution,
% truncated over the region [l,u], where -a<l<u<a for some
% 'a' and l and u are column vectors;
% uses acceptance rejection and inverse-transform method;
% -----------------------------------------------------------------
function x=tn(l,u)
    tol=2; % controls switch between methods
    % threshold can be tuned for maximum speed for each platform
    % case: abs(u-l)>tol, uses accept-reject from randn
    I=abs(u-l)>tol; x=l;
    if any(I)
        tl=l(I); tu=u(I); x(I)=trnd(tl,tu);
    end
    % case: abs(u-l)<tol, uses inverse-transform
    I=~I;
    if any(I)
        tl=l(I); tu=u(I); pl=erfc(tl/sqrt(2))/2; pu=erfc(tu/sqrt(2))/2;
        x(I)=sqrt(2)*erfcinv(2*(pl-(pl-pu).*rand(size(tl))));
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% trnd - uses acceptance rejection to simulate from truncated normal
% -----------------------------------------------------------------
function  x=trnd(l,u)
    x=randn(size(l)); % sample normal
    % keep list of rejected
    I=find(x<l|x>u); d=length(I);
    while d>0 % while there are rejections
        ly=l(I); % find the thresholds of rejected
        uy=u(I);
        y=randn(size(ly));
        idx=y>ly&y<uy; % accepted
        x(I(idx))=y(idx); % store the accepted
        I=I(~idx); % remove accepted from list
        d=length(I); % number of rejected
    end
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% ConvMetric - robust convergence metric function
% -----------------------------------------------------------------
function [err_small,err_wrms,ewt] = ConvMetric(xnew,xold,atol,rtol)
    
    % number of variables
    Nvars = length(xnew);

    % error weights
    ewt = 1./(atol + 0.5*abs(xnew+xold)*rtol);
    
    % error weighted-rms norm
    err_wrms = norm((xnew-xold).*ewt)/sqrt(Nvars);
    %err_wrms = norm(xnew-xold)/(norm(xnew+xold)/2);
    
    % convergence indicator
    err_small = err_wrms <= 1.0;
    %err_small = err_wrms <= rtol;
end
% -----------------------------------------------------------------

% -----------------------------------------------------------------
% PrintScreen - print on screen function
% -----------------------------------------------------------------
function  MyString = PrintScreen(t,F,err,x)

    % number of variables
    Nvars = length(x);

    % print header in the first level
    if t == 1 && Nvars <= 5
        MyHeader = '\n level  func value   error        design variables \n';
        fprintf(MyHeader);
    elseif t == 1 && Nvars > 5
        MyHeader = 'It is not possible to print more than 5 design variables on the screen \n';
        fprintf(MyHeader);
        MyHeader = '\n level  func value   error \n';
        fprintf(MyHeader);
    end
    
    % initial string with (t, F, err)
    MyString = '\n %5g %+6.9f %+6.9f';
    
	% print values on screen
    if Nvars <= 5
        % string with (t, F, err, x)
        for i=1:Nvars
            MyString = strcat(MyString,' %+6.9f');
        end
        % values with x
        fprintf(MyString,t,F,err,x);
    else
        % values without x
        fprintf(MyString,t,F,err);
    end
end
% -----------------------------------------------------------------
