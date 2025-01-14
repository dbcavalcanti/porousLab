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
        function Sl = saturationDegree(~, pc, porousMedia)
            if (pc <= 0.0)
                Sl = 1.0;
            else
                Sl = 1.0 - 1.9722e-11 * pc^(2.4279);
            end
            Sl = max(Sl, porousMedia.Slr);
        end
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function dSldpc = derivativeSaturationDegree(~, pc, ~)
            if (pc <= 0.0)
                dSldpc = 0.0;
            else
                dSldpc = - 2.4279 * 1.9722e-11 * pc^(1.4279);
            end
        end
        
    end
end