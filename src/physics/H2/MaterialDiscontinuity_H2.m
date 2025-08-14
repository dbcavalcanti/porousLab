%% MaterialDiscontinuity_H class
% This class represents a material discontinuity in a porous medium, 
% characterized by its initial aperture and the fluid properties. 
% It provides methods to compute the longitudinal permeability 
% coefficient and compressibility based on the material's properties.
%
%% Methods
% * *longitudinalPermeability*: Computes the longitudinal permeability 
%                               coefficient based on the cubic law.
% * *compressibility*: Computes  the compressibility of the material 
%                      discontinuity based on its aperture and fluid 
%                      properties.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef MaterialDiscontinuity_H2 < Material_H2
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        initialAperture = 0.0;
        leakoff         = 1.0;
        porosity        = 1.0;
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialDiscontinuity_H2(matData)
            this = this@Material_H2(matData)
            if isempty(matData.initialAperture) == false
                this.initialAperture = matData.initialAperture;
            end
            if isempty(matData.leakoff) == false
                this.leakoff = matData.leakoff;
            end
            if isempty(matData.porosity) == false
                this.porosity = matData.porosity;
            end
        end
    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficients based on the
        % cubic's law
        function k = longitudinalPermeability(this)
            w = this.initialAperture();
            k = this.cubicLaw(w);
        end

    end
    %% Static methods
    methods(Static)
        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficient based on the
        % cubic's law
        function kl = cubicLaw(w)
            kl = w*w*w/12.0;
        end
    end
end