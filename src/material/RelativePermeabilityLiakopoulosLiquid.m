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
classdef RelativePermeabilityLiakopoulosLiquid < RelativePermeability  
    %% Properties
    % Parameters taken from OGS-6. 
    % OGS reference:
    % Asadi, R., Ataie-Ashtiani, B. (2015): A Comparison of finite volume
    % formulations and coupling strategies for two-phase flow in deforming
    % porous media. Comput. Geosci., p. 24ff.
    properties (SetAccess = public, GetAccess = public)
        a     = 2.207;
        b     = 1.0121;
        Slmin = 0.2;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityLiakopoulosLiquid()
            this = this@RelativePermeability('liakopoulosLiquid');
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function klr = calculate(this, Sl, porousMedia)
            if (Sl < this.Slmin)
                klr = porousMedia.klrmin;
            elseif (Sl > 1.0)
                klr = 1.0;
            else
                klr = 1.0 - this.a * (1.0 - Sl)^this.b;
                klr = max(klr,porousMedia.klrmin);
            end
        end
        
    end
end