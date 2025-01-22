%% Material_HM class
%
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef Material_HM < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        fluid       = Fluid();
        porousMedia = PorousMedia();
        mechanical  = [];
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_HM(matData)
            this.fluid       = matData.fluid;
            this.porousMedia = matData.porousMedia;
            % Mechanical constitutive behavior
            if strcmp('elastic',matData.porousMedia.mechanical)
                this.mechanical = MechanicalLinearElastic();
            end
        end
    end
    %% Public methods
    methods
        % -----------------------------------------------------------------
        % Evaluate the mechanical constitutive law
        function [stress,D] = mechanicalLaw(this,ip)
            [stress,D] = this.mechanical.eval(this.porousMedia,ip);
        end

        % -----------------------------------------------------------------
        % Get the number of state variables associated with the mechanical
        % constitutive law
        function nstVar = getNumberStateVar(this)
            nstVar = this.mechanical.nstVar;
        end

        % -----------------------------------------------------------------
        % Returns the biot coefficient
        function biot = biotCoeff(this)
            biot = this.porousMedia.biot;
        end

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
            % Compute teh compressibility
            comp = (biot - phi)/Ks + phi/Kf;
        end
    end
end