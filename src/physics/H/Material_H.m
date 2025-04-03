%% Material_H Class
% This class represents a material model that combines fluid and porous 
% media properties. It provides methods to compute material properties 
% such as permeability tensor and compressibility coefficient.
%
%% Methods
% * *permeabilityTensor*: Computes and returns the permeability tensor of 
%                         the material, which is derived from the 
%                         intrinsic permeability matrix of the porous 
%                         media and the fluid viscosity.
% * *compressibilityCoeff*: Computes and returns the compressibility 
%                           coefficient of the material based on Biot's 
%                           coefficient, porosity, solid bulk modulus, 
%                           and fluid bulk modulus.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef Material_H < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        fluid       = Fluid();
        porousMedia = PorousMedia();
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_H(matData)
            this.fluid       = matData.fluid;
            this.porousMedia = matData.porousMedia;
        end
    end
    %% Public methods
    methods

        % -----------------------------------------------------------------
        % Returns the biot coefficient
        function kh = permeabilityTensor(this)
            kh = this.porousMedia.intrinsicPermeabilityMatrix();
            kh = kh / this.fluid.mu;
        end

        % -----------------------------------------------------------------
        % Computes the compressibility coefficient
        function comp = compressibilityCoeff(this)
            % Get material parameters
            biot = this.porousMedia.biot;     % Biot's coefficient
            phi  = this.porousMedia.phi;      % Porosity
            Ks   = this.porousMedia.Ks;       % Solid bulk modulus
            Kf   = this.fluid.K;              % Fluid bulk modulus
            % Compute the compressibility
            comp = (biot - phi)/Ks + phi/Kf;
        end
    end
end