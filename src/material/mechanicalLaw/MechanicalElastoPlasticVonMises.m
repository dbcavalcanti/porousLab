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
        % f = ||s - beta|| - sqrt(2/3) * K(alpha)
        % Radial return
%         function [stress,Dt] = eval(this,material,ip)
% 
%             % Constitutive matrix
%             De = this.elasticConstitutiveMatrix(material,ip);
%             Dt = De;
% 
%             % Trial stress vector
%             stress = De * (ip.strain - ip.plasticstrainOld);
% 
%             % Evaluate the yield condition
%             f = this.yieldCondition(material,ip,stress);
% 
%             % Elastic step
%             if f < 0.0, return, end
% 
%             svmTrial = this.vonMisesStress(stress);
%             Iv = [1.0;1.0;1.0;0.0];
% 
%             % Hardening variables
%             % alphaOld  = this.statevarOld(1);         % Isotropic hardening
%             % betaOld   = this.statevarOld(2:end);     % Backstress
% 
%             % Decompose the stress tensor
%             p  = this.hydrostaticStress(stress);
%             sdtrial = this.deviatoricStress(stress);
% 
%             n = this.flowVector(material,ip,stress);
% 
%             % Relative stress
%             % xi = sd - betaOld;
% 
%             % Initialize variables for the return mapping
%             lambda    = 0.0;
%             iter      = 1;
%             G         = this.shearModulus(material);
% 
%             % Return mapping: closest point projection
%             while (abs(f) > this.returnYieldConditionTol)
% 
%                 % Get the isotropic hardening
%                 H = 0.0;
% 
%                 % Residual derivative
%                 d = -3.0*G - H;
% 
%                 % Update the plastic multiplier
%                 lambda = lambda - f / d;
% 
%                 % Check yield condition
%                 f = svmTrial - 3.0 * G * lambda - material.sy0;
% 
%                 % Update iteration counter
%                 if iter > this.returnMappingMaxIter, break, end
%                 iter = iter + 1;
% 
%             end
% 
%             % Update variables
%             factor = (1.0 - 3.0*G*lambda/svmTrial);
%             sd = sdtrial * factor;
%             stress = p * Iv + sd;
% 
%             % Update the plastic strain
%             ip.plasticstrain = ip.plasticstrainOld + lambda * sqrt(3.0/2.0) * n;
% 
%             % Tangent tensor
%             K = this.bulkModulus(material);
%             Im = Iv * Iv';
%             Id = this.hessianJ2();
%             
%             Dt = K * Im;
%             Dt = Dt + 2.0 * G * (1.0 - 3.0 * G * lambda /svmTrial) * Id;
%             Dt = Dt + 6.0 * G * G * (lambda/svmTrial - 1.0 / (3.0*G)) * (n * n');
% 
%         end

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