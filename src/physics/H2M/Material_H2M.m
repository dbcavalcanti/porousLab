%% Material_H2 class
%
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef Material_H2M < Material_H2    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        mechanical = [];
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_H2M(matData)
            this = this@Material_H2(matData);
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
        % Evaluate the mechanical constitutive law
        function [cul,cug] = mechanicalCompressibilityCoeffs(this,Sl)
            cul = this.porousMedia.biot * Sl;
            cug = this.porousMedia.biot * (1.0 - Sl);
        end

        % -----------------------------------------------------------------
        function flag = hasPlasticStrain(this)
            flag = this.mechanical.isElastoPlastic();
        end
    end
end