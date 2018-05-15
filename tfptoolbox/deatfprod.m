function [ out ] = deatfprod( X, Y, W, P, varargin )
%DEATFPROD Data Envelopment Analysis Price-weighted productivity indices decomposition 
%   This function computes Price-weighted productivity indices decompositions
%   based on Data Envelopment Analysis following  Balk, B.M. and Zofío J.L. (2018). 
%   "The Many Decompositions of Total Factor Productivity Change." 
%   No. ERS-2018-003-LIS. ERIM Report Series Research in Management. 
%   Erasmus Research Institute of Management. Erasmus University, 
%   URL http://hdl.handle.net/1765/104721.
%
%   out = DEATFPROD(X, Y, W, P, Name, Value) computes Data Envelopment Analysis 
%   Price-weighted productivity indices decomposition for inputs X and 
%   outputs Y, with prices W and P respectively. Model properties are specified 
%   using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'orient': orientation. Input orientated 'io' (default), output 
%     orientated 'oo'.
%   - 'period': reference period. Uses geometric mean of two consecutive 
%     periods as reference: 'geomean' (default). To obtain the results 
%     using the base period only, set 'period' to 'base'. To set the 
%     comparison period as reference use 'comparison'.
%   - 'decomp': decomposition. 'complete' (default) following Balk and 
%     Zofío (2018), 'ccd' (Caves, Christensen, and Diewert, 1982), 'crs' 
%     constant returns to scale.
%   - 'names': DMU names.
%
%   Example
%     
%      tfprod_io_geomean_complete = deatfprod(X, Y, W, P, 'orient', 'io',...
%                           'period', 'geomean', 'decomp', 'complete');
%
%   See also DEATFPM, DEATFPMB, DEATFPGPROD
%
%   Copyright 2018 Bert M. Balk, Javier Barbero, Jose L. Zofio
%   http://www.tfptoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 13, May, 2018
%


    % Check size
    if size(X,1) ~= size(Y,1)
        error('Number of rows in X must be equal to number of rows in Y')
    end
    
    if size(X,3) ~= size(Y,3)
        error('Number of time periods in X and Y must be equal')
    end
    
    if size(W) ~= size(X)
        error('Size of W must be equal to size of X')
    end
    
    if size(P) ~= size(Y)
        error('Size of P must be equal to size of Y')
    end
    
    % Get number of DMUs (m), inputs (m) and outputs (s)
    [n, m, T] = size(X);
    s = size(Y,2);
    
    % Default optimization-options
    optimoptsdef = optimoptions('linprog','display','off', 'Algorithm','dual-simplex', 'TolFun', 1e-10, 'TolCon', 1e-7);
  
    % Parse Parameters
    p = inputParser;

    addParameter(p, 'names', cellstr(int2str((1:n)')),...
                @(x) iscellstr(x) && (length(x) == n) );
    addParameter(p, 'optimopts', optimoptsdef, @(x) ~isempty(x));
    addParameter(p, 'orient','oo',...
                @(x) any(validatestring(x,{'io','oo'})));
    addParameter(p, 'period', 'geomean', @(x) any(validatestring(x, {'base','comparison','geomean'})) );
    addParameter(p, 'fixbaset', 0, @(x) ismember(x, [0, 1]));
    addParameter(p, 'decomp', 'complete',...
                @(x) any(validatestring(x, {'complete', 'ccd', 'crs', })));
    addParameter(p, 'warning', 0,  @(x) ismember(x, [0, 1]));
    
    p.parse(varargin{:})
    options = p.Results;
    
    % Set returnsto scale
    switch(options.decomp)
        case {'complete', 'rts', 'ccd'}
            rts = 'vrs';
        case 'crs'
            rts = 'crs';
    end
    
    % Compute VRS or CRS (if CRS decomposition) Malmquist
    Mvrs = deamalm(X, Y, 'orient', options.orient, 'period', options.period, 'fixbaset', options.fixbaset, 'optimopts', options.optimopts,...
            'rts', rts, 'warning', options.warning);
    
    % Set EC and TC to the VRS solution
    EC = Mvrs.eff.MTEC ;
    TC = Mvrs.eff.MTC ;        

    % Generate output structure
    out.n = n;
    out.m = m;
    out.s = s;
    
    out.names = options.names;
    
    out.model = 'tfp-prod';
    out.orient = options.orient;
    out.rts = rts;
    
    out.period = options.period;
    
    % Store common output in eff structure
    if strcmp(rts, 'vrs')
        tfp.Mvrs = Mvrs.eff.M;
    end
    tfp.EC = EC;
    tfp.TC = TC;

    % Create matrices to store results
    SEC     = nan(n, T - 1);
    OME     = nan(n, T - 1);
    IME     = nan(n, T - 1);
    
    SEC_100 = nan(n, T - 1);
    OME_110 = nan(n, T - 1);
    
    SEC_101 = nan(n, T - 1);
    OME_010 = nan(n, T - 1);
    
    SEC_110 = nan(n, T - 1);
    IME_100 = nan(n, T - 1);
    
    SEC_010 = nan(n, T - 1);
    IME_101 = nan(n, T - 1);
    
    Prodcheck = nan(n, T - 1);   
    
    switch(options.period)
        % For base and comparison periods compute the MB index
        case {'base', 'comparison'} 

            % For each time period
            for t=1:T-1
                % Get base period
                if isempty(options.fixbaset) || options.fixbaset == 0
                    tb = t;
                elseif options.fixbaset == 1
                    tb = 1;
                end
                
                % t of the Price-weighted productivity formula
                switch(options.period)
                    case 'base'
                        tmb = tb;
                    case 'comparison'
                        tmb = t + 1;
                end
                
               switch(options.orient)
                   case 'oo'

                       % Scale Effect (Paths A and C)
                       SECd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,tb), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       SECd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       SECd3 = sum(W(:,:,tmb) .* X(:,:,tb), 2);
                       SECd4 = sum(W(:,:,tmb) .* X(:,:,t+1), 2);
                       
                       SECd1.eff = 1 ./ SECd1.eff;
                       SECd2.eff = 1 ./ SECd2.eff;  
      
                        
                       SEC_100(:, t) = (SECd1.eff ./ SECd2.eff) .* (SECd3 ./ SECd4);
                        
                       % Outuput Mix Effect (Paths A and C)
                       OMEd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       OMEd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       OMEd3 = sum(P(:,:,tmb) .* Y(:,:,t+1), 2);
                       OMEd4 = sum(P(:,:,tmb) .* Y(:,:,tb), 2);
                       
                       OMEd1.eff = 1 ./ OMEd1.eff;
                       OMEd2.eff = 1 ./ OMEd2.eff;
                       
                       OME_110(:, t) = (OMEd1.eff ./ OMEd2.eff) .* (OMEd3 ./ OMEd4);
                        
                       % Scale effect (Paths B and D)
                       SECd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       SECd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
             
                       SECd1.eff = 1 ./ SECd1.eff;
                       SECd2.eff = 1 ./ SECd2.eff; 
                       
                       SEC_101(:, t) = (SECd1.eff ./ SECd2.eff) .* (SECd3 ./ SECd4);
                        
                       % Outuput Mix Effect (Paths B and D)
                       OMEd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       OMEd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                        
                       OMEd1.eff = 1 ./ OMEd1.eff;
                       OMEd2.eff = 1 ./ OMEd2.eff;
                       
                       OME_010(:, t) = (OMEd1.eff ./ OMEd2.eff) .* (OMEd3 ./ OMEd4);
                       
                       % Check
                       Prodcheck(:,t) =  (OMEd3 ./ OMEd4) ./ (SECd4 ./ SECd3);
                       
                       % SEC and OME geomean
                       SEC(:,t) = (SEC_100(:, t) .* SEC_101(:, t)).^(1/2);
                       OME(:,t) = (OME_110(:, t) .* OME_010(:, t)).^(1/2);                       
                                               
                   case 'io'
                       % Scale Effect (Paths E and G)
                       SECd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       SECd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       SECd3 = sum(P(:,:,tmb) .* Y(:,:,t+1), 2);
                       SECd4 = sum(P(:,:,tmb) .* Y(:,:,tb), 2);
                       
                       SECd1.eff = 1 ./ SECd1.eff;
                       SECd2.eff = 1 ./ SECd2.eff;
                          
                       SEC_110(:, t) = (SECd1.eff ./ SECd2.eff) .* (SECd3 ./ SECd4);
                        
                       % Input Mix Effect (Paths E and G)
                       IMEd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       IMEd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);                       
                       IMEd3 = sum(W(:,:,tmb) .* X(:,:,tb), 2);
                       IMEd4 = sum(W(:,:,tmb) .* X(:,:,t+1), 2);
                       
                       IMEd1.eff = 1 ./ IMEd1.eff;
                       IMEd2.eff = 1 ./ IMEd2.eff;   
                       
                       IME_100(:, t) = (IMEd1.eff ./ IMEd2.eff) .* (IMEd3 ./ IMEd4);
                       
                       % Scale effect (Paths F and H)
                       SECd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,t+1) , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       SECd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,tb)  , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                        
                       SECd1.eff = 1 ./ SECd1.eff;
                       SECd2.eff = 1 ./ SECd2.eff;
                       
                       SEC_010(:, t) = (SECd1.eff ./ SECd2.eff) .* (SECd3 ./ SECd4);
                       
                       % Inpu Mix Effect (Paths B and D)
                       IMEd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       IMEd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                       
                       IMEd1.eff = 1 ./ IMEd1.eff;
                       IMEd2.eff = 1 ./ IMEd2.eff;
                       
                       IME_101(:, t) = (IMEd1.eff ./ IMEd2.eff) .* (IMEd3 ./ IMEd4);
                       
                       % Check
                       Prodcheck(:,t) =  (SECd3 ./ SECd4 ) ./ ( IMEd4 ./ IMEd3 );
                       
                       % SEC and IME geomean
                       SEC(:,t) = (SEC_110(:, t) .* SEC_010(:, t)).^(1/2);
                       IME(:,t) = (IME_100(:, t) .* IME_101(:, t)).^(1/2);
                    
                end
                              
                  
            end
                        
        case 'geomean'

            % Get Prod for base and comparison periods
            Prodb = deatfprod( X, Y, W, P, varargin{:}, 'period', 'base' ,'warning', options.warning);
            Prodc = deatfprod( X, Y, W, P, varargin{:}, 'period', 'comparison' ,'warning', options.warning);
            
            % Compute geomeans
            switch(options.orient)
                case 'oo'
                    SEC_100 = (Prodb.tfp.SEC_100 .* Prodc.tfp.SEC_100).^(1/2);
                    SEC_101 = (Prodb.tfp.SEC_101 .* Prodc.tfp.SEC_101).^(1/2);
                    
                    OME_110 = (Prodb.tfp.OME_110 .* Prodc.tfp.OME_110).^(1/2);
                    OME_010 = (Prodb.tfp.OME_010 .* Prodc.tfp.OME_010).^(1/2);
                    
                    SEC = (SEC_100 .* SEC_101).^(1/2);
                    OME = (OME_110 .* OME_010).^(1/2);     
                    
                case 'io'
                    SEC_110 = (Prodb.tfp.SEC_110 .* Prodc.tfp.SEC_110).^(1/2);
                    IME_100 = (Prodb.tfp.IME_100 .* Prodc.tfp.IME_100).^(1/2);
                    
                    SEC_010 = (Prodb.tfp.SEC_010 .* Prodc.tfp.SEC_010).^(1/2);
                    IME_101 = (Prodb.tfp.IME_101 .* Prodc.tfp.IME_101).^(1/2);
                    
                    SEC = (SEC_110 .* SEC_010).^(1/2);
                    IME = (IME_100 .* IME_101).^(1/2);
            end
            
            % Prodcheck
            Prodcheck = (Prodb.tfp.Prodcheck .* Prodc.tfp.Prodcheck).^(1./2);
    
        otherwise
            
            error('Period must be ''base'', ''comparison'', or ''geomean''');
           
            
    end
    
    % Price-weighted productivity index
    switch(options.orient)
        case 'oo'
            Prod = EC .* TC .* SEC .* OME;
            
        case 'io'
            Prod = EC .* TC .* SEC .* IME;

    end
    
    % Substitute full expresion computed MB with MBcheck if NaN
    Prod(isnan(Prod)) = Prodcheck(isnan(Prod));
    
    % Productivity index text
    Prodtext = 'Prod = ';
    switch(options.period)
        case 'base'
            Prodtext = [Prodtext, 'Laspeyres'];
        case 'comparison'
            Prodtext = [Prodtext, 'Paasche'];
        case 'geomean'
            Prodtext = [Prodtext, 'Fisher'];
    end
    
    % Display string and output text
    dispstrCommon = 'names/tfp.Prod/';

    switch(options.decomp)
        case {'complete', 'crs'}
            dispstrCommon = [dispstrCommon, '/tfp.EC/tfp.TC'];
        case 'ccd'
            dispstrCommon = [dispstrCommon, '/tfp.Mvrs'];
    end
    
    switch(options.period)
        case {'base','comparison'}
            
            switch(options.orient)
                case 'oo'
                    dispstr = [dispstrCommon, '/tfp.SEC/tfp.OME/tfp.SEC_100/tfp.OME_110/tfp.SEC_101/tfp.OME_010'];
                case 'io'
                    dispstr = [dispstrCommon, '/tfp.IME/tfp.SEC/tfp.IME_100/tfp.SEC_110/tfp.IME_101/tfp.SEC_010'];
            end
        case 'geomean'
            switch(options.orient)
                case 'oo'
                    dispstr = [dispstrCommon, '/tfp.SEC/tfp.OME'];
                case 'io'
                    dispstr = [dispstrCommon, '/tfp.IME/tfp.SEC'];
            end     
               
    end
                
    % Store results in output structure
    tfp.Prod = Prod;   
    
    tfp.SEC = SEC;
    tfp.OME = OME;
    tfp.IME = IME;
    
    tfp.SEC_100 = SEC_100;
    tfp.OME_110 = OME_110;
    
    tfp.SEC_101 = SEC_101;
    tfp.OME_010 = OME_010;
    
    tfp.SEC_110 = SEC_110;
    tfp.IME_100 = IME_100;
    
    tfp.SEC_010 = SEC_010;
    tfp.IME_101 = IME_101;
    
    tfp.Prodcheck = Prodcheck;
        
    % Return eff structure
    out.tfp = tfp;
    
    % Display string
    out.dispstr = dispstr;
    
    % Custom display texts
    out.disptext_title = 'Total Factor Productivity (TFP)';
    out.disptext_text2 = 'Price-weighted productivity index:';
    
    out.disptext_tfp_Prod = 'Prod';
    out.disptext_tfp_EC = 'EC';
    out.disptext_tfp_TC = 'TC';
    out.disptext_tfp_SEC = 'SEC';
    out.disptext_tfp_OME = 'OME';
    out.disptext_tfp_SEC_100 = 'SEC_100';
    out.disptext_tfp_OME_110 = 'OME_110';
    out.disptext_tfp_SEC_101 = 'SEC_101';    
    out.disptext_tfp_OME_010 = 'OME_010';
    out.disptext_tfp_IME = 'IME';
    out.disptext_tfp_IME_100 = 'IME_100';
    out.disptext_tfp_SEC_110 = 'SEC_110';
    out.disptext_tfp_IME_101 = 'IME_101';
    out.disptext_tfp_SEC_010 = 'SEC_010';
    out.disptext_tfp_RTS = 'RTS';
    out.disptext_tfp_Mvrs = 'Mvrs';
    
    switch(options.decomp)
        case 'complete'
            switch(options.orient)
                case 'oo'
                    out.disptext_text4 = [Prodtext, '. EC = Efficiency Change. TC = Technological Change. SEC = Scale Effect. OME = Output Mix Effect.'];
                case 'io'
                    out.disptext_text4 = [Prodtext, '. EC = Efficiency Change. TC = Technological Change. IME = Input Mix Effect. SEC = Scale Effect. '];
            end
        case 'ccd'
            switch(options.orient)
                case 'oo'
                    out.disptext_text4 = [Prodtext, '. Mvrs = Malmquist VRS. SEC = Scale Effect. OME = Output Mix Effect.'];
                case 'io'
                    out.disptext_text4 = [Prodtext, '. Mvrs = Malmquist VRS. IME = Input Mix Effect. SEC = Scale Effect. '];
            end
        case 'crs'
            switch(options.orient)
                case 'oo'
                    out.disptext_text4 = [Prodtext, '. EC = Efficiency Change. TC = Technological Change. SEC = Scale Effect. OME = Output Mix Effect. Under CRS the SEC factor measures only the input mix effect. '];
                case 'io'
                    out.disptext_text4 = [Prodtext, '. EC = Efficiency Change. TC = Technological Change. IME = Input Mix Effect. SEC = Scale Effect. Under CRS the SEC factor measures only the output mix effect. '];
            end            
    end

end

