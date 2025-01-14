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
        function klr = calculate(~, Sl, porousMedia)
            if (Sl <= porousMedia.Slr)
                Sl = porousMedia.Slr;
            elseif (Sl > 1.0)
                Sl = 1.0;
            end
            klr = 1.0 - 2.207 * (1.0 - Sl)^(1.0121);
            klr = max(klr,porousMedia.klrmin);
        end
        
    end
end