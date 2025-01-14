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
classdef RelativePermeabilityBrooksCoreyLiquid < RelativePermeability  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityBrooksCoreyLiquid()
            this = this@RelativePermeability('brooksCoreyLiquid');
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function klr = calculate(~, Sl, porousMedia)
            Se = porousMedia.effectiveSaturationDegree(Sl);
            klr = Se^((2.0 + 3.0 * porousMedia.lambda)/porousMedia.lambda);
            klr = max(klr,porousMedia.klrmin);
        end
        
    end
end