function [ T ] = tfp2table( out, dispstr )
%TFP2TABLE Convert 'tfpout' results into a table object
%   Convert total factor productivity results in a 'tfpout' structure into
%   a MATLAB table object.
%
%   T = TFP2TABLE( out ) Converts 'tfpout' structure into a MATLAB table
%   object.
%   T = TFP2TABLE( out, dispstr ) Converts 'tfpout' structure into a MATLAB
%   table object using the specified 'dispstr' structure.
%
%   Example
%       
%       io = deatfpm(X, Y, 'orient', 'io', 'period', 'geomean');
%       T = tfp2table(io);
%
%       T2 = tfp2table(io, 'names/tfp.M');
%
%   See also TFPDISP, DEATFPM, DEATFPMB, DEATFPPROD, DEATFPGPROD
%
%   Copyright 2018 Bert M. Balk, Javier Barbero, Jose L. Zofio
%   http://www.tfptoolbox.com
%
%   Version: 1.0
%   LAST UPDATE: 22, December, 2017
%

    if nargin < 2
        dispstr = out.dispstr;
    end    
    
    % Convert to table
    dispstr = strsplit(dispstr, '/');
    
    % Create empty table
    T = table();
    
    % Build Header and Body
    for i=1:length(dispstr)  
            % Get param name
            paramstr = char(dispstr(i));
            
            % Get data
            dat = eval(sprintf('out.%s', paramstr));
            
            % Append to Table
            T = [T table(dat)];
            
            % Get variable name
            [name, ~] = getDEAformat(paramstr, out.orient);
            
            % If no name in output structure
            if isempty(name)
                disptext_field = sprintf('disptext_%s', strrep(paramstr,'.','_'));
                if isfield(out, disptext_field)
                    % If custom name exists in the output structure use it
                    name = eval(sprintf('out.%s',disptext_field));
                else
                    % If not, display paramstr name without eff.
                    name = strrep(paramstr, 'eff.', '');
                end
            end
            
            % Store variable name
            T.Properties.VariableNames(size(T,2)) = cellstr(name);
            
    end
    

end

