% -----------------------------------------------------------------
%  ABC.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last updated: April 23, 2022
% -----------------------------------------------------------------
%  This routine employs the ABC algorithm to solve a Bayesian
%  inference problem, which seeks to find the set of parameters
%  that better promotes the agreement between model predictions
%  and a given set of observations.
%  
%  Input:
%  fun    - model function
%  data   - (Ndata x NQoI) reference data
%  lb     - (Nvars x 1) lower bound
%  ub     - (Nvars x 1) upper bound
%  mu     - (Nvars x 1) mean
%  sigma  - (Nvars x 1) standard deviation
%  ABCobj - ABC object struct
%  
%  Output:
%  X_opt  - (Nvars x 1) optimal point
%  F_opt  - optimal value
%  ABCobj - ABC object struct
%  
%  Reference:
%  T. Toni, D. Welch, N. Strelkowa, A. Ipsen, M. Stumpf,
%  Approximate Bayesian computation scheme for parameter 
%  inference and model selection in dynamical systems
%  J. R. Soc. Interface, 6:187-202, 2007
%  https://doi.org/10.1098%2Frsif.2008.0172
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [x_best,y_best,ABCobj] = ABC(fun,data,lb,ub,mu,sigma,ABCobj)

    % check number of arguments
    if nargin < 6
        error('Too few inputs.')
    elseif nargin > 7
        error('Too many inputs.')
    end
    
    % check if lb or ub are empty
    if isempty(lb) || isempty(ub)
        error('lb and ub must be non-empty')
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
    
   % check for consistency in lb and ub
    if size(lb) ~= size(ub)
        error('lb and ub must have the same dimensions')
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
    
	% define mu (if necessary)
    if isempty(mu)
        mu = (ub+lb)/2;
    end
    
    % define sigma (if necessary)
    if isempty(sigma)
        sigma = (ub-lb)/sqrt(12);
    end
    
    % convert mu to a column vector (if necessary)
    [s1,s2] = size(mu);
    if isrow(mu)
        mu = mu(:);
    elseif s2 ~= 1 
        error('mu must be a column vector')
    end
    
    % convert sigma to a column vector (if necessary)
    [s1,s2] = size(sigma);
    if isrow(sigma)
        sigma = sigma(:);
    elseif s2 ~= 1 
        error('sigma must be a column vector')
    end
    
    % check for consistency in mu and sigma
    if size(mu) ~= size(mu)
        error('mu and sigma must have the same dimensions')
    end
    if any(isnan(mu)) || any(~isfinite(mu))
        error('mu must have finite values')
    end
    if any(isnan(sigma)) || any(~isfinite(sigma))
        error('sigma must have finite values')
    end
    if any(mu < lb) || any(mu > ub)
        error('mu must be in [lb,ub] interval')
    end
    if any(sigma <= 0.0)
        error('All components of sigma must be positive')
    end
    
    % number of variables
    Nvars = length(lb);
    
    % check for dimension consistency in x0
    if length(mu) ~= Nvars
        error('x0 must be a Nvars x 1 vector')
    end
    
    % check for dimension consistency in sigma0
    if length(sigma) ~= Nvars
        error('sigma0 must be a Nvars x 1 vector')
    end
    
    % number of data points and outputs
    [Ndata,Nqoi] = size(data);
    
    % model response dimensions
	[Nfun1,Nfun2] = size(fun(zeros(Nvars,1)));
    
    % check for consistency fun dimensions
    if Nfun2 ~= Nqoi
        error('Nfun2 must be equal to Nqoi')
    end
    
    % check if ABCobj is declared
    if nargin == 6
        ABCobj = [];
    end
    
    % number of samples
    if ~isfield(ABCobj, 'Ns')
        ABCobj.Ns= 100;
    elseif ABCobj.Ns < 1
        error('Ns must be a positive integer')
    end
    
    Ns = ABCobj.Ns;

    % tolerance
    if ~isfield(ABCobj, 'tol')
        ABCobj.tol = 1.0e-3;
    elseif ABCobj.tol <= 0.0
        error('tol must be positive')
    end
    
    tol = ABCobj.tol;
    
    % weigths
    if ~isfield(ABCobj, 'weigths')
        ABCobj.weigths = ones(Nqoi,1)/Nqoi;
    elseif size(ABCobj.weigths) ~= [Nqoi,1]
        error('weigths vector dimension must be Nqoi x 1')
    elseif sum(ABCobj.weigths) ~= 1.0
        error('weigths must add one')
    end
    
    weigths = ABCobj.weigths;
    
    % model parameters best value
    x_best = NaN*zeros(Nvars,1);
    
    % model prediction best value
    y_best = NaN*zeros(Nfun1,Nfun2);

    % acceptance/rejection counters
    accept_counter = 0;
    reject_counter = 0;
    
    % minimum value for ABC error function
    J_min = Inf;
    
    % preallocate memory for parameters samples
    x_samples = zeros(Nvars,Ns);
    
    % preallocate memory for accepted/rejected parameters values
    ABCobj.x_accept = zeros(Nvars,Ns);
    ABCobj.x_reject = zeros(Nvars,Ns);
    
    % preallocate memory for accepted predictions values
    ABCobj.y_accept = zeros(Nfun1*Nqoi,Ns);
    ABCobj.y_reject = zeros(Nfun1*Nqoi,Ns);
    
    % limit vectors for standard truncated Gaussian
	supp_l = ((lb - mu)./sigma)*ones(1,Ns);
	supp_u = ((ub - mu)./sigma)*ones(1,Ns);

    % sampling from the prior distribution (truncated Gaussian)
    for ii=1:Nvars
        x_samples(ii,:) = mu(ii) + ...
                          trandn(supp_l(ii,:),supp_u(ii,:))*sigma(ii);
    end
    
    % sampling from the prior distribution (uniform)
	%for ii=1:Nvars
    %    ABCobj.x_samples(ii,:) = unifrnd(lb(ii,1),ub(ii,1),Ns,1);
    %end
    
    % squared norm of data vector
    norm_data_pow2 = data'*data;
    
    % initiate progress bar
    textprogressbar(' ');
    
    % loop for ABC calculation
    for n=1:Ns

        % define model parameters vector
        x = x_samples(:,n);

        % compute the cost function
        y       = fun(x);
        delta_y = data - y(1:Ndata,:);
        J       = sum(weigths.*diag(((delta_y')*delta_y)./norm_data_pow2));
        
        % accept sample if the cost function is lower than a tolerance
        if abs(J) < tol

            % update acceptance counter
            accept_counter = accept_counter + 1;

            % save model parameters accepted values
            ABCobj.x_accept(:,accept_counter) = x;

            % save the model prediction
            ABCobj.y_accept(:,accept_counter) = reshape(y,[Nfun1*Nqoi,1]);

            % check if J is the minimum
            if J < J_min

                % save the best parameters value
                x_best = x;

                % save the best model prediction
                y_best = y;

                % update cost function minimum value
                J_min = J;
            end
        else
            
            % update rejection counter
            reject_counter = reject_counter + 1;
            
            % save model parameters rejected values
            ABCobj.x_reject(:,reject_counter) = x;

            % save the rejected model prediction
            ABCobj.y_reject(:,reject_counter) = reshape(y,[Nfun1*Nqoi,1]);
            
        end
        
        % update progress bar
        textprogressbar((n/Ns)*100);

    end
    
    % space after the progress bar
    textprogressbar(' ');
    fprintf('\n \n');

    % update dimensions of accepted samples matrices
    if accept_counter > 0
        ABCobj.x_accept = ABCobj.x_accept(:,1:accept_counter);
        ABCobj.y_accept = ABCobj.y_accept(:,1:accept_counter);
    end
    
    % update dimensions of rejected samples matrices
    if reject_counter > 0
        ABCobj.x_reject = ABCobj.x_reject(:,1:reject_counter);
        ABCobj.y_reject = ABCobj.y_reject(:,1:reject_counter);
    end
    
    % number of samples
    ABCobj.Ns = Ns;
    
    % acceptance counter
    ABCobj.accept_counter = accept_counter;
    
    % rejection counter
    ABCobj.reject_counter = reject_counter;
    
    % acceptance rate
    ABCobj.accept_rate = accept_counter/Ns;

    % acceptance rate
    ABCobj.reject_rate = reject_counter/Ns;
    
    % error tolerance
    ABCobj.tol = tol;
    
    % minimum value for ABC error function
    ABCobj.J_min = J_min;

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
% trnd - uses acceptance-rejection to simulate from truncated normal
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
function textprogressbar(c)
% This function creates a text progress bar. It should be called 
% with a STRING argument to initialize and terminate. Otherwise, 
% the number corresponding to progress in % should be supplied.
% 
% INPUTS:   C   Either: Text string to initialize or terminate 
%                       Percentage number to show progress 
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m
%
% Author: Paul Proteus 
%         (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version
% 
% Inspired by: 
% http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

    % Initialization
    persistent strCR;           %   Carriage return pesistent variable

    % Vizualization parameters
    strPercentageLength = 10;   %   Length of percentage string (must be >5)
    strDotsMaximum      = 10;   %   The total number of dots in a progress bar

    % Main 

    if isempty(strCR) && ~ischar(c)
        % Progress bar must be initialized with a string
        error('The text progress must be initialized with a string');
    elseif isempty(strCR) && ischar(c)
        % Progress bar - initialization
        fprintf('%s',c);
        strCR = -1;
    elseif ~isempty(strCR) && ischar(c)
        % Progress bar  - termination
        strCR = [];  
        fprintf([c '\n']);
    elseif isnumeric(c)
        % Progress bar - normal progress
        c = floor(c);
        percentageOut = [num2str(c) '%%'];
        percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
        nDots = floor(c/100*strDotsMaximum);
        dotOut = ['[' repmat('=',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
        strOut = [percentageOut dotOut];

        % Print it on the screen
        if strCR == -1
            % Don't do carriage return during first run
            fprintf(strOut);
        else
            % Do it during all the other runs
            fprintf([strCR strOut]);
        end

        % Update carriage return
        strCR = repmat('\b',1,length(strOut)-1);

    else
        % Any other unexpected input
        error('Unsupported argument type');
    end
end
% -----------------------------------------------------------------