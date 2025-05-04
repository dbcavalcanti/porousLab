%% MechanicalElastoPlasticVonMises Class
% This class implements an elasto-plastic constitutive law based on the 
% Von Mises yield criterion. It extends the _MechanicalElastoPlastic_ 
% base class and provides methods for defining the yield condition, 
% flow rule, and hardening behavior.
%
%% Methods
% * *yieldCondition*: Computes the yield condition based on the Von Mises 
%                     stress and the yield stress.
% * *yieldStressGradient*: Computes the gradient of the yield function 
%                          with respect to the stress vector.
% * *flowVector*: Computes the flow direction vector based on the 
%                 deviatoric stress.
% * *flowVectorGradient*: Computes the gradient of the flow vector with 
%                         respect to the stress vector.
% * *hardening*: Returns the hardening modulus from the material 
%                properties.
% * *hardeningStressGradient*: Returns the gradient of the hardening law 
%                              with respect to the stress vector.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
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