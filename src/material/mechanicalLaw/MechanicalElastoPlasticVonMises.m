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
        % Radial return
        function [stress,Dt] = eval(this,material,ip)

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);
            Ce = this.elasticFlexibilityMatrix(material,ip);
            Dt = De;

            % Trial stress vector
            stress = De * (ip.strain - ip.strainOld) + ip.stressOld;

            % Evaluate the yield condition
            f = this.yieldCondition(material,ip,stress);

            % Elastic step
            if f < 0.0, return, end

            % Hardening variables
            alphaOld  = this.statevarOld(1);         % Isotropic hardening
            beta      = this.statevarOld(2:end);     % Backstress

            % Decompose the stress tensor
            p  = this.hydrostaticStress(stress);
            sd = this.deviatoricStress(stress);
            xi = sd - beta;

            % Initialize variables for the return mapping
            lambda    = 0.0;
            ep        = ip.plasticstrainOld;
            epOld     = ip.plasticstrainOld;
            iter      = 1;

            % Return mapping: closest point projection
            while (abs(f) > this.returnYieldConditionTol)

                % Flow vector
                n = xi / norm(xi);

                % Gradient of the yield condition
                df = this.yieldStressGradient(material,ip,stress);

                % Gradient of the flow rule vector
                dn = this.flowVectorGradient(material,ip,stress);
                
                % Hardening
                h = this.hardening(material,ip,stress);

                % Auxiliary matrix
                Psi = Ce + lambda * dn;

                % Increment of the plastic multiplier
                dlambda = (f - df'*(Psi \ r)) / (df' * (Psi \ n) + h);

                % Update the stress vector
                dstress = -Psi \ (r + dlambda * n);
                stress = stress + dstress;

                % Update variables
                lambda = lambda + dlambda;
                ep     = ep - Ce * dstress;
                alpha  = alphaOld + sqrt(2.0/3.0) * lambda;

                % Check yield condition
                f = this.yieldCondition(material,ip,stress);

                % Residual of the flow rule
                r = -ep + epOld + lambda * n;

                % Update iteration counter
                if iter > this.returnMappingMaxIter, break, end
                iter = iter + 1;

            end

            % Update variables
            beta = betaOld + sqrt(2.0/3.0)* lambada;

            % Compute the flow vector at the final stress state
            n  = this.flowVector(material,ip,stress);
            df = this.yieldStressGradient(material,ip,stress);
            dn = this.flowVectorGradient(material,ip,stress);

            % Update the plastic strain
            ip.plasticstrain = ep;

            % Compute algorithmic tangent constitutive tensor
            Psi = inv(Ce + lambda * dn);
            Dt  = Psi - (Psi * (n * df') * Psi)/(df' * (Psi * n));
            % Dt  = De - (De * (n * df') * De)/(df' * (De * n));

        end

        %------------------------------------------------------------------
        % Yield function definition
        function f = yieldCondition(this,material,ip,stress)
            sVM = this.vonMisesStress(stress);
            sy = material.sy0 + material.Kp * ip.statevar;
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
        function h = hardening(~,~,~,~)
            h = 0.0;
        end

        %------------------------------------------------------------------
        % Gradient of the hardening law wrt to the stress vector
        function dh = hardeningStressGradient(~,~,~,~)
            dh = 0.0;
        end

    end
end