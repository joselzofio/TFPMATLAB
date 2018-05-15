%% Data Envelopment Analysis Toolbox examples
%
%   Copyright 2018 Bert M. Balk, Javier Barbero, Jose L. Zofio
%   http://www.tfptoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 14, May, 2018
%

% Clear and clc
clear
clc

%% Section 2: Data structures
% Inputs
X0 = [3 2;4 3;4 5;8 10;6 4];
X1 = [5 3;3 2;4 3;4 8;6 8];
X = X0;
X(:,:,2) = X1;

% Input prices
W0 = [2 1];
W1 = [3 4];
W = W0;
W(:,:,2) = W1;

% Outputs
Y0 = [9 5;12 10;11 10;4 3;6 10];
Y1 = [7 4;11 9;12 10;5 4;4 6];

Y = Y0;
Y(:,:,2) = Y1;

% Output prices
P0 = [3 2];
P1 = [3 5];

P = P0;
P(:,:,2) = P1;

save 'TFPdata'

%% Section 3: Malmquist productivity index
% Section 3.1: The output oriented malmquist productivity index

% The output orientated base period Malmquist Productivity Index.
load 'TFPdata'
tfpm_oo_base_complete = deatfpm(X, Y, 'orient', 'oo',...
                        'period', 'base', 'decomp', 'complete');
tfpdisp(tfpm_oo_base_complete)

% Alternative decompositions of the output orientated Malmquist Productivity Index
tfpm_oo_base_rts = deatfpm(X, Y, 'orient', 'oo',...
                           'period', 'base', 'decomp', 'rts');
tfpdisp(tfpm_oo_base_rts)

tfpm_oo_base_ccd = deatfpm(X, Y, 'orient', 'oo',...
                           'period', 'base', 'decomp', 'ccd');
tfpdisp(tfpm_oo_base_ccd)

rng(1234567); % Set seed for reproducibility
deatestrtsm(X, Y, 'orient', 'oo', 'nreps', 200, 'alpha', 0.05, 'disp', 1);

tfpm_oo_base_crs = deatfpm(X, Y, 'orient', 'oo',...
                           'period', 'base', 'decomp', 'crs');
tfpdisp(tfpm_oo_base_crs)

% The output orientated-comparison period Malmquist Productivity Index.
tfpm_oo_comp_complete = deatfpm(X, Y, 'orient', 'oo',...
                        'period', 'comparison', 'decomp', 'complete');
tfpdisp(tfpm_oo_comp_complete)

% The output orientated-geometric mean Malmquist Productivity Index.
tfpm_oo_geomean_complete = deatfpm(X, Y, 'orient', 'oo',...
                           'period', 'geomean', 'decomp', 'complete');
tfpdisp(tfpm_oo_geomean_complete)

% 3.2. The input orientated Malmquist Productivity Index.

% The input orientated-base period viewpoint
tfpm_io_base_complete = deatfpm(X, Y, 'orient', 'io',...
                        'period', 'base', 'decomp', 'complete');
tfpdisp(tfpm_io_base_complete)

% Alternative decompositions of the input orientated Malmquist Productivity Index.
tfpm_io_base_rts = deatfpm(X, Y, 'orient', 'io',...
                   'period', 'base', 'decomp', 'rts');
tfpdisp(tfpm_io_base_rts)

tfpm_io_base_ccd = deatfpm(X, Y, 'orient', 'io',...
                   'period', 'base', 'decomp', 'ccd');
tfpdisp(tfpm_io_base_ccd)

% The input orientated-comparison period Malmquist Productivity Index.
tfpm_io_comp_complete = deatfpm(X, Y, 'orient', 'io',...
                        'period', 'comparison', 'decomp', 'complete');
tfpdisp(tfpm_io_comp_complete)

% The input orientated-geometric mean Malmquist Productivity Index.
tfpm_io_geomean_complete = deatfpm(X, Y, 'orient', 'io',...
                           'period', 'geomean', 'decomp', 'complete');
tfpdisp(tfpm_io_geomean_complete)

%% 4. Moorsteen-Bjurek Productivity Index.

% 4.1. The output orientated Moorsteen-Bjurek Productivity Index.

% The output orientated-base period decomposition of the Moorsteen-Bjurek PI
tfpmb_oo_base_complete = deatfpmb(X, Y, 'orient', 'oo',...
                         'period', 'base', 'decomp', 'complete');
tfpdisp(tfpmb_oo_base_complete)

% Alternative decompositions of the Moorsteen-Bjurek productivity index
tfpmb_oo_base_ccd = deatfpmb(X, Y, 'orient', 'oo',...
                    'period', 'base', 'decomp', 'ccd');
tfpdisp(tfpmb_oo_base_ccd)

tfpmb_oo_base_complete = deatfpmb(X, Y, 'orient', 'oo',...
                         'period', 'base', 'decomp', 'crs');
tfpdisp(tfpmb_oo_base_complete)

% The output orientated-comparison period decomposition of the Moorsteen-Bjurek PI
tfpmb_oo_comp_complete = deatfpmb(X, Y, 'orient', 'oo',...
                         'period', 'comparison', 'decomp', 'complete');
tfpdisp(tfpmb_oo_comp_complete)

% The output orientated-geometric mean decomposition of the Moorsteen-Bjurek PI
tfpmb_oo_geo_complete = deatfpmb(X, Y, 'orient', 'oo',...
                         'period', 'geomean', 'decomp', 'complete');
tfpdisp(tfpmb_oo_geo_complete)

% 4.2. The input orientated Moorsteen-Bjurek Productivity Index.

% The input orientated-base period decomposition of the Moorsteen-Bjurek PI.
tfpmb_io_base_complete = deatfpmb(X, Y, 'orient', 'io',...
                         'period', 'base', 'decomp', 'complete');
tfpdisp(tfpmb_io_base_complete)

% The input orientated-comparison period decomposition of the Moorsteen-Bjurek PI
tfpmb_io_comp_complete = deatfpmb(X, Y, 'orient', 'io',...
                         'period', 'comparison', 'decomp', 'complete');
tfpdisp(tfpmb_io_comp_complete)

% The input orientated-geometric mean decomposition of the Moorsteen-Bjurek PI
tfpmb_io_geo_complete = deatfpmb(X, Y, 'orient', 'io',...
                         'period', 'geomean', 'decomp', 'complete');
tfpdisp(tfpmb_io_geo_complete)

%% 5. Price-Weighted Productivity Indices (Fisher)

% 5.1 The output orientated price-weighted productivity index.

% The output orientated-base period decomposition of the price-weighted PI}
tfprod_oo_base_complete = deatfprod(X, Y, W, P, 'orient', 'oo',...
                          'period', 'base', 'decomp', 'complete');
tfpdisp(tfprod_oo_base_complete)

% Alternative decompositions of the price weighted  productivity index
tfprod_oo_base_ccd = deatfprod(X, Y, W, P, 'orient', 'oo',...
                     'period', 'base', 'decomp', 'ccd');
tfpdisp(tfprod_oo_base_ccd)

tfprod_oo_base_crs = deatfprod(X, Y, W, P, 'orient', 'oo',...
                     'period', 'base', 'decomp', 'crs');
tfpdisp(tfprod_oo_base_crs)

% The output orientated-comparison period decomposition of the price-weighted PI
tfprod_oo_comp_complete = deatfprod(X, Y, W, P, 'orient', 'oo',...
                          'period', 'comparison', 'decomp', 'complete');
tfpdisp(tfprod_oo_comp_complete)

% The output orientated-geometric mean decomposition of the price-weighted PI
tfprod_oo_geo_complete = deatfprod(X, Y, W, P, 'orient', 'oo',...
                         'period', 'geomean', 'decomp', 'complete');
tfpdisp(tfprod_oo_geo_complete)

% 5.2 The input orientated price-weighted productivity index.

% The input  orientated-base period decomposition of the price-weighted PI
tfprod_io_base_complete = deatfprod(X, Y, W, P, 'orient', 'io',...
                          'period', 'base', 'decomp', 'complete');
tfpdisp(tfprod_io_base_complete)

% The input orientated-comparison period decomposition of the price-weighted PI
tfprod_io_comp_complete = deatfprod(X, Y, W, P, 'orient', 'io',...
                          'period', 'comparison', 'decomp', 'complete');
tfpdisp(tfprod_io_comp_complete)

% The input orientated-geometric mean decomposition of the price-weighted PI
tfprod_io_geo_complete = deatfprod(X, Y, W, P, 'orient', 'io',...
                         'period', 'geomean', 'decomp', 'complete');
tfpdisp(tfprod_io_geo_complete)

%% 6. Share-Weighted Productivity Indices (Törnqvist).

% 6.1 The output orientated share-weighted productivity index.

% The output orientated-base period decomposition of the share-weighted PI.
tfpgprod_oo_base_complete = deatfpgprod(X, Y, W, P, 'orient', 'oo',...
                            'period', 'base', 'decomp', 'complete');
tfpdisp(tfpgprod_oo_base_complete)

% Alternative decompositions of the share-weighted  productivity index
tfpgprod_oo_base_ccd = deatfpgprod(X, Y, W, P, 'orient', 'oo',...
                       'period', 'base', 'decomp', 'ccd');
tfpdisp(tfpgprod_oo_base_ccd)

tfpgprod_oo_base_crs = deatfpgprod(X, Y, W, P, 'orient', 'oo',...
                       'period', 'base', 'decomp', 'crs');
tfpdisp(tfpgprod_oo_base_crs)

% The output orientated-comparison period decomposition of the share-weighted PI
tfpgprod_oo_comp_complete = deatfpgprod(X, Y, W, P, 'orient', 'oo',...
                            'period', 'comparison', 'decomp', 'complete');
tfpdisp(tfpgprod_oo_comp_complete)

% The output orientated-geometric mean decomposition of the share-weighted PI
tfpgprod_oo_geo_complete = deatfpgprod(X, Y, W, P, 'orient', 'oo',...
                           'period', 'geomean', 'decomp', 'complete');
tfpdisp(tfpgprod_oo_geo_complete)

% 6.2 The input orientated share-weighted productivity index.

% The input  orientated-base period decomposition of the share-weighted PI
tfpgprod_io_base_complete = deatfpgprod(X, Y, W, P, 'orient', 'io',...
                            'period', 'base', 'decomp', 'complete');
tfpdisp(tfpgprod_io_base_complete)

% The input orientated-comparison period decomposition of the price-weighted PI
tfpgprod_io_comp_complete = deatfpgprod(X, Y, W, P, 'orient', 'io',...
                            'period', 'comparison', 'decomp', 'complete');
tfpdisp(tfpgprod_io_comp_complete)

% The input orientated-geometric mean decomposition of the share-weighted PI
tfpgprod_io_geo_complete = deatfpgprod(X, Y, W, P, 'orient', 'io',...
                           'period', 'geomean', 'decomp', 'complete');
tfpdisp(tfpgprod_io_geo_complete)

%% 7. Advanced options, displaying and exporting results

% 7.1 Specifying DMU names.
names = {'A'; 'B'; 'C'; 'D'; 'E'};

tfpgprod_names = deatfpgprod(X, Y, W, P, 'orient', 'io',...
                 'period', 'geomean', 'decomp', 'complete',...
                 'names', names);
tfpdisp(tfpgprod_names)

% 7.2 Custom display.
tfpgprod_custom = deatfpgprod(X, Y, W, P, 'orient', 'io',...
                  'period', 'geomean', 'decomp', 'complete');
tfpdisp(tfpgprod_custom, 'names/tfp.GProd')             

% 7.3 Exporting results.
tfpgprod_exp = deatfpgprod(X, Y, W, P, 'orient', 'io',...
                'period', 'geomean', 'decomp', 'complete');
T = tfp2table(tfpgprod_exp);
writetable(T, 'tfpgprod_results.csv');            

                       
                       

