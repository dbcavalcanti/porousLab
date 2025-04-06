%% NonlinearScheme_Picard Class
% This class implements the Picard nonlinear scheme for solving nonlinear 
% systems of equations. It is based on a fully implicit time integration 
% scheme and includes an optional relaxation mechanism to improve 
% convergence.
%
% Reference:
% Li and Wei (2015). An efficient finite element procedure for analyzing 
% three‐phase porous media based on the relaxed Picard method. Int J Numer 
% Methods Eng, 101(11):825-846.
%
%% Methods
% * *assembleLinearSystem*: Assembles the Jacobian matrix and the residual 
%                           vector for the nonlinear system.
% * *applyBCtoRHS*: Applies boundary conditions to the right-hand side 
%                   vector.
% * *addNodalForces*: Adds nodal forces to the right-hand side vector. 
% * *eval*: Evaluates the solution increment and updates the solution 
%           vector.
% * *convergence*: Checks for convergence of the nonlinear solver.
% * *updateRelaxation*: Updates the relaxation parameter based on the 
%                       generalized angle between successive increments.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef NonlinearScheme_Picard < NonlinearScheme
    %% Public properties
    properties(SetAccess = public,GetAccess = public)
        relax = 1.0;
        applyRelaxation = false;
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = NonlinearScheme_Picard()
            this = this@NonlinearScheme();
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Assembles the linear system
        function [A,b] = assembleLinearSystem(~,C,K,~,fe,dfidx,~,xOld,dt)
            % RHS vector
            b = C * xOld / dt + fe;

            % LHS matrix
            A = K + dfidx + C / dt;
        end

        %------------------------------------------------------------------
        % Applies the boundary conditions to the right-hand side vector
        function bf = applyBCtoRHS(~,A,b,x,doffree,doffixed)
            bf = b(doffree) - A(doffree,doffixed) * x(doffixed);
        end

        %------------------------------------------------------------------
        % Adds nodal forces to the right-hand side vector
        function b = addNodalForces(~,b,fe)
            b = b + fe;
        end

        %------------------------------------------------------------------
        % Evaluates the solution increment and updates the solution vector
        function [X,dx] = eval(this,A,b,X,dxOld,freedof,iter)
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
        % Checks the convergence of the nonlinear solver
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
        % Updates the relaxation parameter based on the generalized angle 
        % between successive increments.
        function updateRelaxation(this,X,XOld,dxOld)
            dx = X - XOld;

            % Compute generalized angle between successive increments
            gamma = acos((dx' * dxOld) / norm(dx) * norm(dxOld));

            % Update relaxation parameter
            if (gamma < pi/4.0)
                this.relax = min(this.relax * 2.0, 1.0);
            elseif (gamma > pi/2.0)
                this.relax = this.relax / 2.0;
            end
        end
    end
end
