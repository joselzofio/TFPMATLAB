%% Total Factor Productivity Toolbox examples
%
%   Copyright 2018 Bert M. Balk, Javier Barbero, Jose L. Zofio
%   http://www.tfptoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 3, July, 2018
%

% Clear and clc
clear
clc

%% Section 2: Data structures

% Data for this section
load 'DataAgriculture'

% Inputs
X1960 = [CAPITAL1960, LAND1960, LABOR1960, INTINP1960];
X2004 = [CAPITAL2004, LAND2004, LABOR2004, INTINP2004];

% Outputs
Y1960 = [LS1960, CROPS1960, FARMOUT1960];
Y2004 = [LS2004, CROPS2004, FARMOUT2004];

% Inputs prices
W1960 = [WCAPITAL1960, WLAND1960, WLABOR1960, WINTINP1960];
W2004 = [WCAPITAL2004, WLAND2004, WLABOR2004, WINTINP2004];

% Output prices
P1960 = [PLS1960, PCROPS1960, PFARMOUT1960];
P2004 = [PLS2004, PCROPS2004, PFARMOUT2004];

% Merge the two years in a 3D array
X = X1960;
X(:,:,2) = X2004;

Y = Y1960;
Y(:,:,2) = Y2004;

W = W1960;
W(:,:,2) = W2004;

P = P1960;
P(:,:,2) = P2004;

X = X ./ 1000;
Y = Y ./ 1000;

%% Section 3:  Malmquist Productivity Index

% Section 3.1: The output orientated MPI
tfpm_oo_base_complete = deatfpm(X, Y, 'orient', 'oo', 'period', 'base',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpm_oo_base_complete)        

tfpm_oo_base_rts = deatfpm(X, Y, 'orient', 'oo', 'period', 'base',...
                        'decomp', 'rts', 'names', STATE_NAME);
tfpdisp(tfpm_oo_base_rts)

tfpm_oo_base_ccd = deatfpm(X, Y, 'orient', 'oo', 'period', 'base',...
                        'decomp', 'ccd', 'names', STATE_NAME);
tfpdisp(tfpm_oo_base_ccd)

rng(1234567); % Set seed for reproducibility
deatestrtsm(X, Y, 'orient', 'oo', 'nreps', 200, 'alpha', 0.05, 'disp', 1);

tfpm_oo_base_crs = deatfpm(X, Y, 'orient', 'oo', 'period', 'base',...
                        'decomp', 'crs', 'names', STATE_NAME);
tfpdisp(tfpm_oo_base_crs)

tfpm_oo_comp_complete = deatfpm(X, Y, 'orient', 'oo', 'period', 'comparison',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpm_oo_comp_complete),

tfpm_oo_geomean_complete = deatfpm(X, Y, 'orient', 'oo', 'period', 'geomean',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpm_oo_geomean_complete)

% Section 3.2: The input orientated MPI
tfpm_io_base_complete = deatfpm(X, Y, 'orient', 'io', 'period', 'base',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpm_io_base_complete)

%% Section 4: Moorsteen-Bjurek Productivity Index

% Section 4.1: The output orientated decomposition of MBPI
tfpmb_oo_base_complete = deatfpmb(X, Y, 'orient', 'oo', 'period', 'base',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpmb_oo_base_complete)

tfpmb_oo_geo_complete = deatfpmb(X, Y, 'orient', 'oo', 'period', 'geomean',...
                        'decomp', 'complete','names', STATE_NAME);
tfpdisp(tfpmb_oo_geo_complete)

% Section 4.2: The input orientated decomposition of MBPI
tfpmb_io_base_complete = deatfpmb(X, Y, 'orient', 'io', 'period', 'base',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpmb_io_base_complete)

%% Section 5: Price-Weighted Productivity Indices (Fisher)

% Section 5.1: The output orientated decomposition of a price-weighted productivity index
tfprod_oo_base_complete = deatfprod(X, Y, W, P, 'orient', 'oo', 'period', 'base',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfprod_oo_base_complete)

tfprod_oo_comp_complete = deatfprod(X, Y, W, P, 'orient', 'oo', 'period', 'comparison',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfprod_oo_comp_complete)

tfprod_oo_geo_complete = deatfprod(X, Y, W, P, 'orient', 'oo', 'period', 'geomean',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfprod_oo_geo_complete)

% Section 5.2: The input orientated decomposition of a price-weighted productivity index
tfprod_io_base_complete = deatfprod(X, Y, W, P, 'orient', 'io', 'period', 'base',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfprod_io_base_complete)

%% Section 6: Share-Weighted Productivity Indices (Tornqvist)

% Section 6.1: The output orientated decomposition of a share-weighted productivity index
tfpgprod_oo_base_complete = deatfpgprod(X, Y, W, P, 'orient', 'oo', 'period', 'base',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpgprod_oo_base_complete)

tfpgprod_oo_comp_complete = deatfpgprod(X, Y, W, P, 'orient', 'oo', 'period', 'comparison',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpgprod_oo_comp_complete)

tfpgprod_oo_geo_complete = deatfpgprod(X, Y, W, P, 'orient', 'oo', 'period', 'geomean',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpgprod_oo_geo_complete)

% Section 6.2: The input orientated decomposition of a share-weighted productivity index
tfpgprod_io_base_complete = deatfpgprod(X, Y, W, P, 'orient', 'io', 'period', 'base',...
                        'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpgprod_io_base_complete)

%% Section 7: Displaying and exporting results

% Section 7.1: Custom display
tfpgprod_custom = deatfpgprod(X, Y, W, P, 'orient', 'io', 'period', 'geomean', ...
                  'decomp', 'complete', 'names', STATE_NAME);
tfpdisp(tfpgprod_custom, 'names/tfp.GProd')

% Section 7.2: Exporting results
tfpgprod_exp = deatfpgprod(X, Y, W, P, 'orient', 'io', 'period', 'geomean',...
               'decomp', 'complete', 'names', STATE_NAME);
T = tfp2table(tfpgprod_exp);
writetable(T, 'tfpgprod_results.csv');

