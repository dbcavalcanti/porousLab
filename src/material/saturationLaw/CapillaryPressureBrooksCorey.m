%% CapillaryPressureBrooksCorey class
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
classdef CapillaryPressureBrooksCorey < CapillaryPressure  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = CapillaryPressureBrooksCorey()
            this = this@CapillaryPressure('brooksCorey');
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function Sl = saturationDegree(~, pc, porousMedia)
            if (pc >= porousMedia.Pb)
                Se = (porousMedia.Pb/pc)^porousMedia.lambda;
                Sl = Se * (1.0 - porousMedia.Slr - porousMedia.Sgr) + porousMedia.Slr;
            else
                Sl = 1.0;
            end
            Sl = max(min(Sl, 1.0 - porousMedia.Sgr - eps), porousMedia.Slr + eps);
        end
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function dSldpc = derivativeSaturationDegree(this, pc, porousMedia)
            % if (pc >= porousMedia.Pb)
            % 
            % else
            %     dSldpc = 1.0;
            % end
            Sl = this.saturationDegree(pc, porousMedia);
            dPcdSl = (porousMedia.Pb / (porousMedia.lambda * (porousMedia.Slr - Sl))) * ((Sl - porousMedia.Slr)/(1.0 - porousMedia.Sgr - porousMedia.Slr)) ^ (-1.0/porousMedia.lambda);
            dSldpc = 1.0 / dPcdSl;
        end
        
    end
end