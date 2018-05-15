function [ ] = tfpdisp(out, dispstr)
%TFPDISP Display total factor productivity results
%   Display total factor productivity results stored in a 'tfpout'
%   structure.
%   TFPDISP( out ) Display total factor productivity results.
%   TFPDISP( out, dispstr ) Display results using the specified 'dispstr'.
%
%   Example
%       
%      tfpm_io_geomean_complete = deatfpm(X, Y, 'orient', 'io',...
%                           'period', 'geomean', 'decomp', 'complete');
%
%      tfpdisp(tfpm_io_geomean_complete, 'names/tfp.M');
%
%   See also TFP2TABLE
%
%   Copyright 2018 Bert M. Balk, Javier Barbero, Jose L. Zofio
%   http://www.tfptoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 14, May, 2018
%

    if nargin < 2
        deadisp(out);
    else
        deadisp( out, dispstr );
    end
        
end

