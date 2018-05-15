function [ out ] = deatfpmb( X, Y, varargin )
%DEATFPMB Data Envelopment Analysis Moorsteen-Bjurek indices decomposition 
%   This function computes Moorsteen-Bjurek indices decompositions based on
%   Data Envelopment Analysis following  Balk, B.M. and Zofío J.L. (2018). 
%   "The Many Decompositions of Total Factor Productivity Change." 
%   No. ERS-2018-003-LIS. ERIM Report Series Research in Management. 
%   Erasmus Research Institute of Management. Erasmus University, 
%   URL http://hdl.handle.net/1765/104721.
%
%   out = DEATFPMB(X, Y, Name, Value) computes Data Envelopment Analysis 
%   Moorsteen-Bjurek indices decomposition for inputs X and outputs Y. Model 
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
%     Zofío (2018), 'ccd' (Caves, Christensen, and Diewert, 1982), 'crs' 
%     constant returns to scale.
%   - 'names': DMU names.
%
%   Example
%     
%      tfpmb_io_geomean_complete = deatfpmb(X, Y, 'orient', 'io',...
%                           'period', 'geomean', 'decomp', 'complete');
%
%   See also DEATFP, DEATFPROD, DEATFPGPROD
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
    
    out.model = 'tfp-mb';
    out.orient = options.orient;
    out.rts = rts;
    
    out.period = options.period;
    
    % Store common output in eff structure
    tfp.Mvrs = Mvrs.eff.M;
    tfp.EC = EC;
    tfp.TC = TC;
    
    % Decompositions
    switch(options.decomp)
        case {'complete', 'ccd', 'crs'}
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

        MBcheck = nan(n, T - 1);


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

                    % Compute mean values
                    Xbar = (X(:,:,tb) + X(:,:,t + 1)) ./ 2 ;
                    Ybar = (Y(:,:,tb) + Y(:,:,t + 1)) ./ 2 ;

                   switch(options.orient)
                       case 'oo'

                           % Scale Effect (Paths A and C)
                           SECd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,tb), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           SECd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           SECd3 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Ybar     , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           SECd4 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Ybar     , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);

                           SECd1.eff = 1 ./ SECd1.eff;
                           SECd2.eff = 1 ./ SECd2.eff;  
                           SECd3.eff = 1 ./ SECd3.eff; 
                           SECd4.eff = 1 ./ SECd4.eff; 

                           SEC_100(:, t) = (SECd1.eff ./ SECd2.eff) .* (SECd3.eff ./ SECd4.eff);

                           % Outuput Mix Effect (Paths A and C)
                           OMEd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', Xbar      , 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           OMEd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', Xbar      , 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           OMEd3 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           OMEd4 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);

                           OMEd1.eff = 1 ./ OMEd1.eff;
                           OMEd2.eff = 1 ./ OMEd2.eff;
                           OMEd3.eff = 1 ./ OMEd3.eff;
                           OMEd4.eff = 1 ./ OMEd4.eff;

                           OME_110(:, t) = (OMEd1.eff ./ OMEd2.eff) .* (OMEd3.eff ./ OMEd4.eff);

                           % Scale effect (Paths B and D)
                           SECd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           SECd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);

                           SECd1.eff = 1 ./ SECd1.eff;
                           SECd2.eff = 1 ./ SECd2.eff; 

                           SEC_101(:, t) = (SECd1.eff ./ SECd2.eff) .* (SECd3.eff ./ SECd4.eff);

                           % Outuput Mix Effect (Paths B and D)
                           OMEd3 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           OMEd4 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);

                           OMEd3.eff = 1 ./ OMEd3.eff;
                           OMEd4.eff = 1 ./ OMEd4.eff;

                           OME_010(:, t) = (OMEd1.eff ./ OMEd2.eff) .* (OMEd3.eff ./ OMEd4.eff);

                           % Check
                           MBcheck(:,t) =  (OMEd1.eff ./ SECd4.eff) ./ (OMEd2.eff ./ SECd3.eff);

                           % SEC and OME geomean
                           SEC(:,t) = (SEC_100(:, t) .* SEC_101(:, t)).^(1/2);
                           OME(:,t) = (OME_110(:, t) .* OME_010(:, t)).^(1/2);     

                       case 'io'

                           % Scale Effect (Paths E and G)
                           SECd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           SECd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           SECd3 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', Xbar      , 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           SECd4 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', Xbar      , 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'oo', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);

                           SECd1.eff = 1 ./ SECd1.eff;
                           SECd2.eff = 1 ./ SECd2.eff;
                           SECd3.eff = 1 ./ SECd3.eff;
                           SECd4.eff = 1 ./ SECd4.eff;

                           SEC_110(:, t) = (SECd1.eff ./ SECd2.eff) .* (SECd3.eff ./ SECd4.eff);

                           % Input Mix Effect (Paths E and G)
                           IMEd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Ybar      , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           IMEd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Ybar      , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           IMEd3 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           IMEd4 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,tb) , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);

                           IMEd1.eff = 1 ./ IMEd1.eff;
                           IMEd2.eff = 1 ./ IMEd2.eff;
                           IMEd3.eff = 1 ./ IMEd3.eff;
                           IMEd4.eff = 1 ./ IMEd4.eff;

                           IME_100(:, t) = (IMEd1.eff ./ IMEd2.eff) .* (IMEd3.eff ./ IMEd4.eff);

                           % Scale effect (Paths F and H)
                           SECd1 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,t+1) , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           SECd2 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb), 'Yeval', Y(:,:,tb)  , 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);

                           SECd1.eff = 1 ./ SECd1.eff;
                           SECd2.eff = 1 ./ SECd2.eff;

                           SEC_010(:, t) = (SECd1.eff ./ SECd2.eff) .* (SECd3.eff ./ SECd4.eff);

                           % Inpu Mix Effect (Paths B and D)
                           IMEd3 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,t+1), 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);
                           IMEd4 = dea(X(:,:,tmb), Y(:,:,tmb), 'Xeval', X(:,:,tb) , 'Yeval', Y(:,:,t+1), 'rts', rts, 'orient', 'io', 'optimopts', options.optimopts, 'secondstep', 0 ,'warning', options.warning);

                           IMEd3.eff = 1 ./ IMEd3.eff;
                           IMEd4.eff = 1 ./ IMEd4.eff;

                           IME_101(:, t) = (IMEd1.eff ./ IMEd2.eff) .* (IMEd3.eff ./ IMEd4.eff);

                           % Check
                           MBcheck(:,t) =  (SECd3.eff ./ IMEd2.eff ) ./ ( SECd4.eff ./ IMEd1.eff  );

                           % SEC and IME geomean
                           SEC(:,t) = (SEC_110(:, t) .* SEC_010(:, t)).^(1/2);
                           IME(:,t) = (IME_100(:, t) .* IME_101(:, t)).^(1/2);

                   end

                end

            case 'geomean'

                % Get MB for base and comparison periods
                MBb = deatfpmb( X, Y, varargin{:}, 'period', 'base' ,'warning', options.warning);
                MBc = deatfpmb( X, Y, varargin{:}, 'period', 'comparison' ,'warning', options.warning);

                % Compute geomeans
                switch(options.orient)
                    case 'oo'
                        SEC_100 = (MBb.tfp.SEC_100 .* MBc.tfp.SEC_100).^(1/2);
                        SEC_101 = (MBb.tfp.SEC_101 .* MBc.tfp.SEC_101).^(1/2);

                        OME_110 = (MBb.tfp.OME_110 .* MBc.tfp.OME_110).^(1/2);
                        OME_010 = (MBb.tfp.OME_010 .* MBc.tfp.OME_010).^(1/2);

                        SEC = (SEC_100 .* SEC_101).^(1/2);
                        OME = (OME_110 .* OME_010).^(1/2);  

                    case 'io'
                        SEC_110 = (MBb.tfp.SEC_110 .* MBc.tfp.SEC_110).^(1/2);
                        IME_100 = (MBb.tfp.IME_100 .* MBc.tfp.IME_100).^(1/2);

                        SEC_010 = (MBb.tfp.SEC_010 .* MBc.tfp.SEC_010).^(1/2);
                        IME_101 = (MBb.tfp.IME_101 .* MBc.tfp.IME_101).^(1/2);

                        SEC = (SEC_110 .* SEC_010).^(1/2);
                        IME = (IME_100 .* IME_101).^(1/2);
                end
                
                % Prodcheck
                MBcheck = (MBb.tfp.MBcheck .* MBc.tfp.MBcheck).^(1./2);

            otherwise

                error('Period must be ''base'', ''comparison'', or ''geomean''');


        end

        % Moorsteen-Bjureck index
        switch(options.orient)
            case 'oo'
                MB = EC .* TC .* SEC .* OME;

            case 'io'
                MB = EC .* TC .* SEC .* IME;

        end

        % Substitute full expresion computed MB with MBcheck if NaN
        MB(isnan(MB)) = MBcheck(isnan(MB));

        % Display string and output text
        dispstrCommon = 'names/tfp.MB/';

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
        
    end
    
    % Store results in output structure
    tfp.MB = MB;  
    
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
    
    tfp.MBcheck = MBcheck;
        
    % Return eff structure
    out.tfp = tfp;
    
    % Display string
    out.dispstr = dispstr;
    
    % Custom display texts
    out.disptext_title = 'Total Factor Productivity (TFP)';
    out.disptext_text2 = 'Moorsteen-Bjurek index:';
    
    out.disptext_tfp_MB = 'MB';
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
                    out.disptext_text4 = 'MB = Moorsteen-Bjurek. EC = Efficiency Change. TC = Technological Change. SEC = Scale Effect. OME = Output Mix Effect.';
                case 'io'
                    out.disptext_text4 = 'MB = Moorsteen-Bjurek. EC = Efficiency Change. TC = Technological Change. IME = Input Mix Effect. SEC = Scale Effect. ';
            end
        case 'ccd'
            switch(options.orient)
                case 'oo'
                    out.disptext_text4 = 'MB = Moorsteen-Bjurek. Mvrs = Malmquist VRS. SEC = Scale Effect. OME = Output Mix Effect.';
                case 'io'
                    out.disptext_text4 = 'MB = Moorsteen-Bjurek. Mvrs = Malmquist VRS. IME = Input Mix Effect. SEC = Scale Effect. ';
            end
        case 'crs'
            switch(options.orient)
                case 'oo'
                    out.disptext_text4 = 'MB = Moorsteen-Bjurek. EC = Efficiency Change. TC = Technological Change. SEC = Scale Effect. OME = Output Mix Effect. Under CRS the SEC factor measures only the input mix effect. ';
                case 'io'
                    out.disptext_text4 = 'MB = Moorsteen-Bjurek. EC = Efficiency Change. TC = Technological Change. IME = Input Mix Effect. SEC = Scale Effect. Under CRS the SEC factor measures only the output mix effect. ';
            end            
    end

end

