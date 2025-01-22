%% Brooks and Corey for the gas phase
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef RelativePermeabilityPolynomialGas < RelativePermeability  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityPolynomialGas()
            this = this@RelativePermeability('polynomialGas');
        end
    end

    %% Public methods
    methods
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function kgr = calculate(~, Sl, porousMedia)
            kgr = (1.0 - Sl)*porousMedia.m;
            kgr = max(kgr,porousMedia.kgrmin);
        end
        
    end
end