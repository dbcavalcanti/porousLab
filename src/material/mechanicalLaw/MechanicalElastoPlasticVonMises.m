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
classdef MechanicalElastoPlasticVonMises < MechanicalElastoPlastic  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalElastoPlasticVonMises()
            this = this@MechanicalElastoPlastic();
            this.nstVar = 5;   % Hardening + Kinematic hardening
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Yield function definition
        function f = yieldCondition(this,material,ip,stress)
            sVM = this.vonMisesStress(stress);
            sy = material.sy0;% + material.Kp * ip.statevar(1);
            f  = sVM - sy;
        end

        %------------------------------------------------------------------
        % Gradient of the yield function wrt to the stress vector
        function df = yieldStressGradient(this,~,~,stress)
            df = this.vonMisesStressGradient(stress);
        end

        %------------------------------------------------------------------
        % Flow vector
        function n = flowVector(this,~,~,stress)
            % n = this.vonMisesStressGradient(stress);
            sd = this.deviatoricStress(stress);
            n = sd / norm(sd);
        end

        %------------------------------------------------------------------
        % Flow vector gradient
        function dn = flowVectorGradient(this,~,~,stress)
            dn = this.vonMisesStressHessian(stress);
        end

        %------------------------------------------------------------------
        % Hardening law
        function h = hardening(~,material,~,~)
            h = material.Kp;
        end

        %------------------------------------------------------------------
        % Gradient of the hardening law wrt to the stress vector
        function dh = hardeningStressGradient(~,~,~,~)
            dh = 0.0;
        end

    end
end