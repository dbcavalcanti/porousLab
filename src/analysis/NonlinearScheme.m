%% NonlinearScheme class
%
% This in an abstract class that defines a NonlinearScheme object.
%
classdef NonlinearScheme < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        tol = 1.0e-5;
        normalizeError = false;
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = NonlinearScheme()
            return;
        end
    end

    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
        [A,b] = assembleLinearSystem(this,C,K,fi,fe,dfidx,x,xOld,dt);

        %------------------------------------------------------------------
        bf = applyBCtoRHS(this,A,b,x,doffree,doffixed);

        %------------------------------------------------------------------
        b = addNodalForces(this,b,fe);

        %------------------------------------------------------------------
        [X,dx] = eval(this,J,r,X,dx,freedof,iter);

        %------------------------------------------------------------------
        convFlg = convergence(this,X,XOld,dx,b,doffree,iter);
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        function setConvergenceTolerance(this,tol)
            this.tol = tol;
        end
    end
end
