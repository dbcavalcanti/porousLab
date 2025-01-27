%% MechanicalElastoPlastic class
%
% This class defines a generic elastoplastic constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef MechanicalElastoPlastic < MechanicalLinearElastic  
    properties (SetAccess = public, GetAccess = public)
        returnMappingMaxIter = 100;
        returnYieldConditionTol = 1.0e-8;
        returnNormResidualFlowRuleTol = 1.0e-8;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalElastoPlastic()
            this = this@MechanicalLinearElastic();
        end
    end
    %% Abstract methods
    methods(Abstract)

        f = yieldCondition(this,material,ip,stress);

        df = yieldStressGradient(this,material,ip,stress);

        n = flowVector(this,material,ip,stress);

        dn = flowVectorGradient(this,material,ip,stress);

        h = hardening(this,material,ip,stress);

        dh = hardeningStressGradient(this,material,ip,stress);

    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = eval(this,material,ip)

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);
            Ce = this.elasticFlexibilityMatrix(material,ip);
            Dt = De;

            % Trial stress vector
            stress = De * (ip.strain - ip.plasticstrain);

            % Evaluate the yield condition
            f = this.yieldCondition(material,ip,stress);

            % Elastic step
            if f < 0.0, return, end

            % Initialize variables for the return mapping
            lambda = 0.0;
            ep     = ip.plasticstrainOld;
            epOld  = ip.plasticstrainOld;
            iter   = 1;
            r      = zeros(4,1);

            % Return mapping: closest point projection
            while (abs(f) > this.returnYieldConditionTol) || (norm(r) > this.returnNormResidualFlowRuleTol)

                % Flow vector
                n = this.flowVector(material,ip,stress);

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

                % Update the plastic multipler
                lambda = lambda + dlambda;

                % Update the plastic strain
                ep = ep - Ce * dstress;

                % Check yield condition
                f = this.yieldCondition(material,ip,stress);

                % Residual of the flow rule
                r = -ep + epOld + lambda * n;

                % Update iteration counter
                if iter > this.returnMappingMaxIter, break, end
                iter = iter + 1;

            end

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
    end
    methods (Static)
        function flag = isElastoPlastic()
            flag = true;
        end
    end
end