function [ S, SB, pvalue, critval ] = deatestrtsm( X, Y, varargin )
%DEATESTRTS Data envelopment analysis test of returns to scale
%   Computes data envelopment analysis test of returns to scale (RTS) for
%   multiple time periods.
%
%   [S, SB, pvalue, critval] = DEATESTRTSM(X, Y, Name, Value) computes data 
%   envelopment analysis test of returns to scale with inputs X and outputs
%   Y. Model properties are specified using one or more Name ,Value pair 
%   arguments. The function returns the test statistic (S), the bootstrapped
%   statistic (SB), the p-value (pvalue) and the critical value (critval).
%
%   Additional properties:
%   - 'orient': orientation. Input oriented 'io', output oriented 'oo'.
%   - 'names': DMU names.
%   - 'nreps': number of bootstrap replications. Default is 200.
%   - 'alpha': alpha value for confidence intervals. Default is 0.05.
%   - 'disp': set to 1 to display results on screen. Default is 0.
%
%   Example
%     
%      [S, SB, pvalue, critval] = deatestrts(X, Y, 'orient', 'io');
%
%   See also DEAOUT, DEA, DEASCALE, DEABOOT
%
%   Copyright 2018 Bert M. Balk, Javier Barbero, Jose L. Zofio
%   http://www.tfptoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 14, January, 2018
%
        
    % Check size
    if size(X,1) ~= size(Y,1)
        error('Number of rows in X must be equal to number of rows in Y')
    end
    
    if size(X,3) ~= size(Y,3)
        error('Number of time periods in X and Y must be equal')
    end
    

    % Parse Parameters
    p = inputParser;
    addParameter(p, 'orient','io',...
                @(x) any(validatestring(x,{'io','oo'})));
    addParameter(p, 'nreps', 200, @(x) isnumeric(x) & (x > 0) );
    addParameter(p, 'alpha', 0.05, @(x) isnumeric(x) & (x > 0));
    addParameter(p, 'disp', 0, @(x) ismember(x, [0, 1]));    
    p.parse(varargin{:})
    options = p.Results;
    
    % Get options
    nreps = options.nreps;
    alph = options.alpha;
    
    % For each time period
    t = size(X,3);
    
    S = nan(t,1);
    SB = nan(t,nreps);
    pvalue = nan(t,1);
    critval = nan(t,1);
    
    for i=1:t
        [St, SBt, pvaluet, critvalt] = deatestrts(X(:,:,i), Y(:,:,i), varargin{:}, 'disp', 0);
                
        S(i,1) = St;
        SB(i,1:nreps) = SBt;
        pvalue(i,1) = pvaluet;
        critval(i,1) = critvalt;
    end
    
    % Display test results
    if options.disp
        fprintf('_______________________________\n');
        fprintf('<strong>DEA Test of RTS</strong>\n\n');
        
        fprintf('H0: Globally CRS \n')
        fprintf('H1: VRS \n\n')
        
        fprintf('Bootstrap replications: %i \n', nreps);
        fprintf('Significance level: %4.2f \n \n', alph);
        
        for i=1:t
            fprintf('<strong>Period %i</strong> \n', i)
            fprintf('S statistic: %7.4f \n', S(i));
            fprintf('Critical value: %7.4f \n', critval(i));
            fprintf('p-value: %7.4f \n', pvalue(i))
        end
        
    end
       

end

