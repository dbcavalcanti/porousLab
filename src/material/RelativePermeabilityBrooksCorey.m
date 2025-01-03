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
classdef RelativePermeabilityBrooksCorey < RelativePermeability  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityBrooksCorey()
            this = this@RelativePermeability('brooksCorey');
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function klr = liquidRelativePermeability(~, Sl, porousMedia)
            Se = porousMedia.effectiveSaturationDegree(Sl);
            klr = Se^((2.0 + 3.0 * porousMedia.lambda)/porousMedia.lambda);
            klr = max(klr,porousMedia.klrmin);
        end
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function kgr = gasRelativePermeability(~, Sl, porousMedia)
            Se = porousMedia.effectiveSaturationDegree(Sl);
            kgr = (1.0 - Se)*(1.0 - Se)*(1.0 - Se^((2.0 + porousMedia.lambda)/porousMedia.lambda));
            kgr = max(kgr,porousMedia.kgrmin);
        end
        
    end
end