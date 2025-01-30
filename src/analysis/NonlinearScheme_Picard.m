%% NonlinearScheme_Newton Class
%
% Considers a fully implicit time integration scheme.
%
% Reference:
% Li, W., & Wei, C. (2015). An efficient finite element procedure for
% analyzing three‐phase porous media based on the relaxed Picard method.
% International Journal for Numerical Methods in Engineering, 101(11),
% 825-846.
%
classdef NonlinearScheme_Picard < NonlinearScheme 
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = NonlinearScheme_Picard()
            this = this@NonlinearScheme();
        end
    end
    %% Public methods
    % Implementation of the abstract methods declared in super-class
    methods (Static)
        %------------------------------------------------------------------
        function [X, dx] = eval(A,b,X,freedof)
            XOld = X;
            X(freedof) = A\b;
            dx = X - XOld;
        end

        %------------------------------------------------------------------
        function [A, b] = assembleLinearSystem(C, K, ~, fe, dfidx, ~, xOld, dt)
            % RHS vector
            b = C * xOld / dt + fe;
            % LHS matrix
            A = K + dfidx + C / dt;
        end

        %------------------------------------------------------------------
        function bf = applyBCtoRHS(A, b, x, doffree, doffixed)
            bf = b(doffree) - A(doffree,doffixed)*x(doffixed);
        end

        %------------------------------------------------------------------
        function b = addNodalForces(b,fe)
            b = b + fe;
        end

    end

    methods

        %------------------------------------------------------------------
        function convFlg = convergence(this,X,~,dX,~, doffree,iter)
            % Evaluate error
            normError = norm(dX(doffree));
            if this.normalizeError
                normError = normError / norm(X(doffree));
                fprintf("\t\t iter.: %3d , ||dX||/||X|| = %7.3e \n",iter,normError);
            else
                fprintf("\t\t iter.: %3d , ||dX||/||X|| = %7.3e \n",iter,normError);
            end
            % Check convergence
            if (norm(normError) < this.tol) && (iter > 1)
                convFlg = true;
            else
                convFlg = false;
            end
        end

    end
    
end