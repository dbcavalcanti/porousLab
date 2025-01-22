%% Brooks and Corey for the gas phase
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
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