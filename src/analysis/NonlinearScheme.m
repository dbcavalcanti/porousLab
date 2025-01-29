%% NonlinearScheme class
%
% This in an abstract class that defines the NonlinearScheme object
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: January 2025
%%%
%
classdef NonlinearScheme < handle
    %% Properties
    properties (SetAccess = private, GetAccess = public)
        tol = 1.0e-5;
    end
    %% Constructor method
    methods
        function this = NonlinearScheme()
        end
    end
    %% Public methods
    methods
        function setConvergenceTolerance(this,tol)
            this.tol = tol;
        end
    end
    %% Abstract methods
    methods (Abstract)
        [X, dx] = eval(J,r,X,freedof);
        [A,b] = assembleLinearSystem(C, K, fi, fe, dfidx, x, xOld, dt);
        bf = applyBCtoRHS(A, b, X, doffree, doffixed);
        convFlg = convergence(X,XOld,dx,b,iter);
    end
end