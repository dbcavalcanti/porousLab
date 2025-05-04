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
classdef MaterialDiscontinuity_H < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        initialAperture = 0.0;
        fluid           = Fluid();
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialDiscontinuity_H(matData)
            this.initialAperture = matData.initialAperture;
            this.fluid = matData.fluid;
        end
    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficient based on the
        % cubic's law
        function kl = longitudinalPermeability(this)
            w = this.initialAperture();
            kl = w*w*w/12.0/this.fluid.mu;
        end

        %------------------------------------------------------------------
        % Computes  the compressibility of the material discontinuity 
        % based on its aperture and fluid properties.
        function c = compressibility(this)
            w = this.initialAperture();
            c = w / this.fluid.K;
        end


    end
end