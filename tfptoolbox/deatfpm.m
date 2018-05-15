function [ out ] = deatfpm( X, Y, varargin )
%DEATFPM. Data Envelopment Analysis Malmquist indices decomposition 
%   This function computes Malmquist indices decompositions based on Data 
%   Envelopment Analysis following  Balk, B.M. and Zofío J.L. (2018). 
%   "The Many Decompositions of Total Factor Productivity Change." 
%   No. ERS-2018-003-LIS. ERIM Report Series Research in Management. 
%   Erasmus Research Institute of Management. Erasmus University, 
%   URL http://hdl.handle.net/1765/104721.
%
%   out = DEATFPM(X, Y, Name, Value) computes Data Envelopment Analysis 
%   Malmquist indices decomposition for inputs X and outputs Y. Model 
%   properties are specified using one or more Name ,Value pair arguments.
%
%   Additional properties:
%   - 'orient': orientation. Input orientated 'io' (default), output 
%     orientated 'oo'.
%   - 'period': reference period. Uses geometric mean of two consecutive 
%     periods as reference: 'geomean' (default). To obtain the results 
%     using the base period only, set 'period' to 'base'. To set the 
%     comparison period as reference use 'comparison'.
%   - 'decomp': decomposition. 'complete' (default) following Balk and 
%     Zofío (2018), 'rts' returns to scale (Ray and Desli, 1997), 'ccd' 
%     (Caves, Christensen, and Diewert, 1982), 'crs' constant returns to 
%     scale.
%   - 'names': DMU names.
%
%   Example
%     
%      tfpm_io_geomean_complete = deatfpm(X, Y, 'orient', 'io',...
%                           'period', 'geomean', 'decomp', 'complete');
%
%   See also DEATFPMB, DEATFPROD, DEATFPGPROD
%
%   Copyright 2018 Bert M. Balk, Javier Barbero, Jose L. Zofio
%   http://www.tfptoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 19, March, 2018
%

    % Check size
    if size(X,1) ~= size(Y,1)
        error('Number of rows in X must be equal to number of rows in Y')
    end
    
    if size(X,3) ~= size(Y,3)
        error('Number of time periods in X and Y must be equal')
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
    addParameter(p, 'orient','io',...
                @(x) any(validatestring(x,{'io','oo'})));
    addParameter(p, 'period', 'geomean', @(x) any(validatestring(x, {'base','comparison','geomean'})) );
    addParameter(p, 'fixbaset', 0, @(x) ismember(x, [0, 1]));
    addParameter(p, 'decomp', 'complete',...
                @(x) any(validatestring(x, {'complete', 'rts', 'ccd', 'crs', })));
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
    
    % Compute CRS and VRS Malmquist
    Mcrs = deamalm(X, Y, 'orient', options.orient, 'period', options.period, 'fixbaset', options.fixbaset, 'optimopts', options.optimopts,...
            'rts', 'crs', 'warning', options.warning);
    Mvrs = deamalm(X, Y, 'orient', options.orient, 'period', options.period, 'fixbaset', options.fixbaset, 'optimopts', options.optimopts,...
            'rts', 'vrs', 'warning', options.warning);
    
    % Set M to the CRS solution
    M = Mcrs.eff.M ;
    
    % Set MTEC and MTC to the VRS solution
    EC = Mvrs.eff.MTEC ;
    TC = Mvrs.eff.MTC ;        
        
    % Compute returns to sacle
    RTS = M ./ Mvrs.eff.M ;

    % Generate output structure
    out.n = n;
    out.m = m;
    out.s = s;
    
    out.names = options.names;
    
    out.model = 'tfp-m';
    out.orient = options.orient;
    out.rts = rts;
    
    out.period = options.period;
    
    % Store M, MTEC, MTC and MRTS in eff structure
    tfp.M = M;
    tfp.Mvrs = Mvrs.eff.M;
    tfp.EC = EC;
    tfp.TC = TC;
    tfp.RTS = RTS;
    
    % Decompositions
    switch(options.decomp)
        case {'complete', 'ccd'}
            % Balk (2001) decomposition
            
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
                        
                        % t of the Moorsteen-Bjureck formula
                        switch(options.period)
                            case 'base'
                                tmb = tb;
                            case 'comparison'
                                tmb = t + 1;
                        end
                         
                       switch(options.orient)
                           case 'oo'
        
                               % Scale Effect (Paths A and C)
                               OSE_10 = deascale(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb), 'orient', 'oo', 'optimopts', options.optimopts, 'warning', options.warning) ;
                               OSE_00 = deascale(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,tb), 'orient', 'oo', 'optimopts', options.optimopts, 'warning', options.warning) ;
                                
                               SEC_100(:, t) = (OSE_10.eff.scale ./ OSE_00.eff.scale);
                                
                               % Outuput Mix Effect (Paths A and C)
                               OSE_11 = deascale(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'orient', 'oo', 'optimopts', options.optimopts, 'warning', options.warning) ;
                       
                               OME_110(:, t) = (OSE_11.eff.scale ./ OSE_10.eff.scale);
                                
                               % Scale effect (Paths B and D)
                               OSE_01 = deascale(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,t+1), 'orient', 'oo', 'optimopts', options.optimopts, 'warning', options.warning) ;
                               
                               SEC_101(:, t) = (OSE_11.eff.scale ./ OSE_01.eff.scale);
                                
                               % Outuput Mix Effect (Paths B and D)
                               OME_010(:, t) = (OSE_01.eff.scale ./ OSE_00.eff.scale);
                               
                                % Inver SEC and OME if orientation is 'oo' (as scale efficiencies with 'oo' are > 1)
                               SEC_100(:, t) = 1 ./ SEC_100(:, t) ;
                               OME_110(:, t) = 1 ./ OME_110(:, t) ;
                               SEC_101(:, t) = 1 ./ SEC_101(:, t) ;
                               OME_010(:, t) = 1 ./ OME_010(:, t) ;
 
                               % SEC and OME geomean
                               SEC(:,t) = (SEC_100(:, t) .* SEC_101(:, t)).^(1/2);
                               OME(:,t) = (OME_110(:, t) .* OME_010(:, t)).^(1/2);     
                                
                           case 'io'

                               % Scale Effect (Paths E and G)
                               ISE_11 = deascale(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'orient', 'io', 'optimopts', options.optimopts, 'warning', options.warning) ;
                               ISE_10 = deascale(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1) , 'Yeval', Y(:,:,tb), 'orient', 'io', 'optimopts', options.optimopts, 'warning', options.warning) ;
                                                              
                               SEC_110(:, t) = (ISE_11.eff.scale ./ ISE_10.eff.scale);
                                
                               % Input Mix Effect (Paths E and G)
                               ISE_10 = deascale(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb), 'orient', 'io', 'optimopts', options.optimopts, 'warning', options.warning) ;
                               ISE_00 = deascale(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,tb), 'orient', 'io', 'optimopts', options.optimopts, 'warning', options.warning) ;
                                  
                               IME_100(:, t) = (ISE_10.eff.scale ./ ISE_00.eff.scale);
                               
                               % Scale effect (Paths F and H)
                               ISE_01 = deascale(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,t+1), 'orient', 'io', 'optimopts', options.optimopts, 'warning', options.warning) ;
                               
                               SEC_010(:, t) = (ISE_01.eff.scale ./ ISE_00.eff.scale);
                               
                               % Inpu Mix Effect (Paths B and D)                               
                               IME_101(:, t) = (ISE_11.eff.scale ./ ISE_01.eff.scale);
                                                       
                               % SEC and IME geomean
                               SEC(:,t) = (SEC_110(:, t) .* SEC_010(:, t)).^(1/2);
                               IME(:,t) = (IME_100(:, t) .* IME_101(:, t)).^(1/2);
                            
                        end
                        
                                                
                    end
                                
                case 'geomean'
                    
                    % Get MB for base and comparison periods
                    Mb = deatfpm( X, Y, varargin{:}, 'period', 'base' ,'warning', options.warning);
                    Mc = deatfpm( X, Y, varargin{:}, 'period', 'comparison' ,'warning', options.warning);
                    
                    % Compute geomeans
                    switch(options.orient)
                        case 'oo'
                            SEC_100 = (Mb.tfp.SEC_100 .* Mc.tfp.SEC_100).^(1/2);
                            SEC_101 = (Mb.tfp.SEC_101 .* Mc.tfp.SEC_101).^(1/2);
                            
                            OME_110 = (Mb.tfp.OME_110 .* Mc.tfp.OME_110).^(1/2);
                            OME_010 = (Mb.tfp.OME_010 .* Mc.tfp.OME_010).^(1/2);
                            
                            SEC = (SEC_100 .* SEC_101).^(1/2);
                            OME = (OME_110 .* OME_010).^(1/2);  
                            
                        case 'io'
                            SEC_110 = (Mb.tfp.SEC_110 .* Mc.tfp.SEC_110).^(1/2);
                            IME_100 = (Mb.tfp.IME_100 .* Mc.tfp.IME_100).^(1/2);
                            
                            SEC_010 = (Mb.tfp.SEC_010 .* Mc.tfp.SEC_010).^(1/2);
                            IME_101 = (Mb.tfp.IME_101 .* Mc.tfp.IME_101).^(1/2);
                            
                            SEC = (SEC_110 .* SEC_010).^(1/2);
                            IME = (IME_100 .* IME_101).^(1/2);
                    end
            
                otherwise
                    
                    error('Period must be ''base'', ''comparison'', or ''geomean''');
                   
                    
            end

            
            % Display string and output text
            dispstrCommon = 'names/tfp.M/';
            
            switch(options.decomp)
                case 'complete'
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
                
        case 'rts'
                        
            % Display string
            dispstr = 'names/tfp.M/tfp.EC/tfp.TC/tfp.RTS';
            
        case 'crs'
            
            % Substitute EC and TC vrs for crs
            tfp.EC = Mcrs.eff.MTEC;
            tfp.TC = Mcrs.eff.MTC;
            
            % Display string
            dispstr = 'names/tfp.M/tfp.EC/tfp.TC';           
                        

    end
    
    % Return eff structure
    out.tfp = tfp;
    
    % Display string
    out.dispstr = dispstr;
    
    % Custom display texts
    out.disptext_title = 'Total Factor Productivity (TFP)';
    out.disptext_text2 = 'Malmquist:';
    
    out.disptext_tfp_M = 'M';
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
                    out.disptext_text4 = 'M = Malmquist. EC = Efficiency Change. TC = Technological Change. SEC = Scale Effect. OME = Output Mix Effect.';
                case 'io'
                    out.disptext_text4 = 'M = Malmquist. EC = Efficiency Change. TC = Technological Change. IME = Input Mix Effect. SEC = Scale Effect. ';
            end
        case 'rts'
            out.disptext_text4 = 'M = Malmquist. EC = Efficiency Change. TC = Technological Change. RTS = Returns to Scale.';
        case 'ccd'
            switch(options.orient)
                case 'oo'
                    out.disptext_text4 = 'M = Malmquist. Mvrs = Malmquist VRS. SEC = Scale Effect. OME = Output Mix Effect.';
                case 'io'
                    out.disptext_text4 = 'M = Malmquist. Mvrs = Malmquist VRS. IME = Input Mix Effect. SEC = Scale Effect. ';
            end
        case 'crs'
            out.disptext_text4 = 'M = Malmquist. EC = Efficiency Change. TC = Technological Change.';
    end
    

end

