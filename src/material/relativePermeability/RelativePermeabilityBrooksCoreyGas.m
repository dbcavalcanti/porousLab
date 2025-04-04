%% RelativePermeabilityBrooksCoreyGas Class
% This class implements the Brooks-Corey model for calculating the 
% relative permeability of the gas phase in porous media. The model is 
% based on the effective saturation degree and incorporates a minimum 
% relative permeability threshold.
%
%% Method
% * *calculate*: Computes the gas phase relative permeability based on the  
%                effective saturation degree Se and the Brooks-Corey model.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef RelativePermeabilityBrooksCoreyGas < RelativePermeability  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityBrooksCoreyGas()
            this = this@RelativePermeability('brooksCoreyGas');
        end
    end

    %% Public methods
    methods
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function kgr = calculate(~, Sl, porousMedia)
            Se = porousMedia.effectiveSaturationDegree(Sl);
            kgr = (1.0 - Se)*(1.0 - Se)*(1.0 - Se^((2.0 + porousMedia.lambda)/porousMedia.lambda));
            kgr = max(kgr,porousMedia.kgrmin);
        end
        
    end
end