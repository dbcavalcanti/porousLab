%% Material_Elastic class
%
% This class defines an linear elastic stress-strain constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef CapillaryPressureLiakopoulos < CapillaryPressure  
    %% Properties
    % Parameters taken from OGS-6. 
    % OGS reference:
    % Asadi, R., Ataie-Ashtiani, B. (2015): A Comparison of finite volume
    % formulations and coupling strategies for two-phase flow in deforming
    % porous media. Comput. Geosci., p. 24ff.
    properties (SetAccess = public, GetAccess = public)
        a     = 1.9722e-11;
        b     = 2.4279;
        Slmin = 0.2;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = CapillaryPressureLiakopoulos()
            this = this@CapillaryPressure('liakopoulos');
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function Sl = saturationDegree(this, pc, ~)
            if (pc < 0.0)
                Sl = 1.0;
            else
                Sl = 1.0 - this.a * pc^(this.b);
            end
            Sl = max(Sl, this.Slmin);
        end
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function dSldpc = derivativeSaturationDegree(this, pc, ~)
            if (pc < 0.0)
                dSldpc = 0.0;
            else
                pcmax  = ((1.0 - this.Slmin)/this.a)^(1.0/this.b);
                pc     = min(pc,pcmax);
                dSldpc = - this.a * this.b * pc^(this.b - 1.0);
            end
        end
        
    end
end