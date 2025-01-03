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
                Sl = porousMedia.Slr + (1.0 - porousMedia.Slr) * (porousMedia.Pb/pc)^porousMedia.lambda;
            else
                Sl = 1.0;
            end
        end
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function dSldpc = derivativeSaturationDegree(~, pc, porousMedia)
            if (pc >= porousMedia.Pb)
                dSldpc = -(1.0 - porousMedia.Slr) * porousMedia.lambda * ((porousMedia.Pb/pc)^porousMedia.lambda) / pc;
            else
                dSldpc = 1.0;
            end
        end
        
    end
end