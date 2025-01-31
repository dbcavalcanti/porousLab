%% NonlinearScheme_Picard Class
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
    properties(SetAccess = public,GetAccess = public)
        relax = 1.0;
        applyRelaxation = true;
    end
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
        function [X, dx] = eval(this,A,b,X,dxOld,freedof,iter)
            XOld = X;
            X(freedof) = A\b;
            if this.applyRelaxation
                if iter > 1
                    this.updateRelaxation(X,XOld,dxOld);
                    X(freedof) = this.relax * X(freedof) + (1.0 - this.relax) * XOld(freedof);
                end
            end
            dx = X - XOld;
        end

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

        %------------------------------------------------------------------
        function updateRelaxation(this,X,XOld,dxOld)
            dx = X - XOld;
            % Compute the generalized angle between successive increments
            gamma = acos((dx'*dxOld)/norm(dx)*norm(dxOld));
            % Update relaxation parameter
            if (gamma < pi/4.0)
                this.relax = min(this.relax * 2.0,1.0);
            elseif (gamma > pi/2.0)
                this.relax = this.relax / 2.0;
            end
        end

    end
    
end