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
classdef MechanicalElastoPlastic < MechanicalLinearElastic  
    properties (SetAccess = public, GetAccess = public)
        returnMappingMaxIter = 20;
        returnYieldConditionTol = 1.0e-8;
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
            Dt = De;

            % Trial stress vector
            stress = De * (ip.strain - ip.strainOld) + ip.stressOld;

            % Evaluate the yield condition
            f = this.yieldCondition(material,ip,stress);

            % Elastic step
            if f < 0.0, return, end

            % Initialize variables for the return mapping
            lambda = 0.0;
            strial = stress;
            iter   = 1;

            % Return mapping
            while f > this.returnYieldConditionTol

                % Flow vector
                n = this.flowVector(material,ip,stress);

                % Gradient of the yield condition
                df = this.yieldStressGradient(material,ip,stress);

                % Gradient of the flow rule vector
                dn = this.flowVectorGradient(material,ip,stress);
                
                % Residual of the flow rule
                r = stress - strial * lambda * De * n;

                % Hardening
                h = this.hardening(material,ip,stress);

                % Auxiliary matrix
                Psi = eye(4) + lambda * De * dn;

                % Increment of the plastic multiplier
                dlambda = (f - df'*(Psi \ r)) / (df' * (Psi \ (De * n)) + h);

                % Update the stress vector
                dstress = -Psi \ (r + dlambda * De * n);
                stress = stress + dstress;

                % Update the plastic multipler
                lambda = lambda + dlambda;

                % Check yield condition
                f = this.yieldCondition(material,ip,stress);

                % Update iteration counter
                if iter > this.returnMappingMaxIter, break, end
                iter = iter + 1;

            end

            % Compute the flow vector at the final stress state
            n  = this.flowVector(material,ip,stress);
            dn = this.flowVectorGradient(material,ip,stress);
            h  = this.hardening(material,ip,stress);

            % Update the plastic strain
            ip.plasticStrain = ip.plasticStrainOld + lambda * n;

            % Compute algorithmic tangent constitutive tensor
            Psi = eye(4) + lambda * De * dn;
            H   = Psi \ De;
            Dt  = H - (H * (n * df') * H)/(df' * (H * n) + h);

        end
    end
    methods (Static)
        function flag = isElastoPlastic()
            flag = true;
        end
    end
end