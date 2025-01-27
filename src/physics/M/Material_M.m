%% Material_M class
%
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef Material_M < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        porousMedia = PorousMedia();
        mechanical  = [];
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_M(matData,lc)
            this.porousMedia = matData.porousMedia;
            % Mechanical constitutive behavior
            if strcmp('elastic',matData.porousMedia.mechanical)
                this.mechanical = MechanicalLinearElastic();
            elseif strcmp('vonMises',matData.porousMedia.mechanical)
                this.mechanical = MechanicalElastoPlasticVonMises();
            elseif strcmp('isoDamage',matData.porousMedia.mechanical)
                this.mechanical = MechanicalIsotropicDamage(lc);
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

        function flag = hasPlasticStrain(this)
            flag = this.mechanical.isElastoPlastic();
        end

    end
end