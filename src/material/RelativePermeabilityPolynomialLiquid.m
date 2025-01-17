%% Brooks and Corey for the gas phase
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef RelativePermeabilityPolynomialLiquid < RelativePermeability  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityPolynomialLiquid()
            this = this@RelativePermeability('polynomialLiquid');
        end
    end

    %% Public methods
    methods
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function klr = calculate(~, Sl, porousMedia)
            klr = Sl*porousMedia.m;
            klr = max(klr,porousMedia.klrmin);
        end
        
    end
end