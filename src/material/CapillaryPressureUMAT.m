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
classdef CapillaryPressureUMAT < CapillaryPressure 
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        pc_curve = [];
        Sl_curve = [];
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = CapillaryPressureUMAT(pc_curve,Sl_curve)
            this = this@CapillaryPressure('umat');
            this.pc_curve = pc_curve;
            this.Sl_curve = Sl_curve;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function Sl = saturationDegree(this, pc, ~)
            if (pc > this.pc_curve(1))
                Sl = this.Sl_curve(1);
            elseif  (pc < this.pc_curve(end))
                Sl = this.Sl_curve(end);
            else
                Sl = interp1(this.pc_curve,this.Sl_curve,pc);
            end
        end
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function dSldpc = derivativeSaturationDegree(this, pc, porousMedia)
            % Pertubation for the numerical derivative
            h = 0.0001;
            % Compute the saturation degree at the perturbed values
            pc_back = this.saturationDegree(pc-h, porousMedia);
            pc_forw = this.saturationDegree(pc+h, porousMedia);
            dSldpc = (pc_forw - pc_back)/(2.0*h);
        end
        
    end
end