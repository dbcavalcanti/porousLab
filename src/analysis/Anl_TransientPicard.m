%% Anl_TransientPicard Class
%
% Performs a transient analysis combined with the Picard method to solve
% the nonlinear system.
%
% Reference:
% Li, W., & Wei, C. (2015). An efficient finite element procedure for
% analyzing three‐phase porous media based on the relaxed Picard method.
% International Journal for Numerical Methods in Engineering, 101(11),
% 825-846.
%
classdef Anl_TransientPicard < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        theta = 1.0;    % Implicit time integration scheme
        dt    = 0.001;  % Time increment;
        tf    = 1.0;    % Final time
        tinit = 0.01;
        dtMax = 0.001;
        dtMin = 0.001;
        adaptStep = false;
        maxIter = 250;
        maxAttempts = 10;
        applyRelax = false;
        relax = 1.0;
        useRelativeError = false;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_TransientPicard(result)
            this = this@Anl('TransientPicard',result);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function process(this,mdl)

            % Initialize transient analysis parameters
            t        = this.tinit;
            t0       = this.tinit;
            maxSteps = length(this.tinit:this.dt:this.tf);
            step     = 1;

            % Initialize solution vector
            XOld = mdl.U;
            X0 = mdl.U;
            X    = XOld;

            % Initialize the output vectors
            this.result.time = zeros(maxSteps,1);
            this.result.p    = zeros(maxSteps,1);
            attempt = 1;
            brokenStep = false;

            % Transient analysis
            while (t0 < this.tf)
                fprintf("\t Time: %12.5f s \n",t);

                % Update transient solution
                X0(mdl.doffree)   = X(mdl.doffree);
                XOld(mdl.doffree) = X(mdl.doffree);

                % Iterative solution
                convFlg = false;
                attemptOld = attempt;
                attempt    = 1;
                
                while attempt < this.maxAttempts
                    fprintf("\t Attempt %3d: \n",attempt);
                    iter = 1;
                    % Iterative process
                    while true
    
                        % Compute model global matrices
                        [K, C, ~, Fext] = mdl.globalMatrices(X);
        
                        % Set transient system
                        A =  K + C / this.dt;
        
                        % Right-handside vector
                        b = Fext + C * X0 / this.dt;

                        % % Apply BC
                        b(mdl.doffree) = b(mdl.doffree) - A(mdl.doffree,mdl.doffixed)*X(mdl.doffixed);

                        % Solve linear system
                        X(mdl.doffree) = A(mdl.doffree,mdl.doffree)\b(mdl.doffree);
                        if (this.applyRelax)
                            if (iter > 1)
                                this.computePicardRelaxation(X, XOld, DX);
                            end
                            X(mdl.doffree) = (1.0 - this.relax)*XOld(mdl.doffree) + this.relax*X(mdl.doffree);
                        end
    
                        % Update variables
                        DX = X - XOld;
                        XOld(mdl.doffree) = X(mdl.doffree);

                        % Compute the residual norm
                        [absError, errorDX] = this.evaluateError(DX,X,mdl);

                        % Print result
                        fprintf("\t\t iter.: %3d , Dp = %7.3e , Dpg = %7.3e , ||Rdp|| = %7.3e \n",iter,absError(1),absError(2),errorDX);
    
                        % Check convergence
                        if (errorDX < 1.0e-5) && (iter > 1)
                            convFlg = true; 
                            break
                        end
    
                        % Check maximum number of iterations
                        iter = iter + 1;
                        if (iter > this.maxIter)
                            break
                        end
                    end

                    % Check convergence
                    if convFlg == true, break, end

                    % Reduce the time step 
                    this.dt = max(this.dt / 4.0, this.dtMin);
                    
                    % Clean the previous attempt
                    XOld  = X0;

                    % Update counter variable
                    attempt = attempt + 1;
                end

                if convFlg == false
                    disp("Transient analysis: solution did not converge!");
                    break
                end

                % Update the state variables
                mdl.updateStateVar();

                % Update time step
                if (this.adaptStep == true) && (attempt == 1) && (brokenStep == false) && (attemptOld == 1)
                    this.dt = min(this.dt * 2,this.dtMax);
                end

                % Save results
                this.result.time(step) = t;
                if isempty(this.result.coordP) == false
                    elemPlot = findElementInMesh(mdl.NODE, mdl.ELEM,this.result.coordP);
                    gle = mdl.element(elemPlot).type.gle;
                    this.result.p(step) = mdl.element(elemPlot).type.pressureField(this.result.coordP,X(gle));
                else
                    this.result.p(step) = 0.0;
                end  
                
                % Update time
                t0 = t;
                if (t + this.dt) > this.tf
                    this.dt = this.tf - t;
                end
                t = t + this.dt;
                step = step + 1;

            end
            
            mdl.updateStateVar();

            % Save final result
            mdl.U = X;
            for el = 1:mdl.nelem
                gle = mdl.element(el).type.gle;
                mdl.element(el).type.ue = mdl.U(gle);
            end

            % Remove unused steps in the output vectors
            this.result.p = this.result.p(1:(step-1));
            this.result.time = this.result.time(1:(step-1));

        end

        %------------------------------------------------------------------
        function [absError, normError] = evaluateError(this,DX,X,mdl)

            absError = [max(abs(DX(mdl.pFreeDof))) , max(abs(DX(mdl.pgFreeDof)))];
            normError = norm(DX(mdl.doffree));
            if (this.useRelativeError)
                normError = normError / norm(X(mdl.doffree));
            end

        end

        %------------------------------------------------------------------
        function printStep(~,X,mdl)

            for i = 1:mdl.nnodes
                fprintf("  %4d: \t",i);
                for j = 1:mdl.ndof_nd
                    fprintf("  %8.4f ",X(mdl.ID(i,j)))
                end
                fprintf("\n");
            end

        end

        %------------------------------------------------------------------
        function [normRp , normRpg] = decomposeResidualVct(~,r, Fext,mdl)

            rp = r(mdl.pFreeDof);
            rpg = r(mdl.pgFreeDof);

            nfep = norm(Fext(mdl.pFreeDof));
            nfepg = norm(Fext(mdl.pgFreeDof));

            if nfep < 1.0e-15, nfep = 1.0; end
            if nfepg < 1.0e-15, nfep = 1.0; end

            normRp = norm(rp)/nfep;
            normRpg = norm(rpg)/nfepg; 

        end

        %------------------------------------------------------------------
        function setUpTransientSolver(this,tinit,dt,tf,dtMax,dtMin,adaptStep)
            if nargin == 4
                dtMax = dt;
                dtMin = dt;
                adaptStep = false;
            end
            this.tinit = tinit;
            this.dt = dt;
            this.tf = tf;
            this.adaptStep = adaptStep;
            this.dtMax = dtMax;
            this.dtMin = dtMin;
        end

        %------------------------------------------------------------------
        function setPicardRelaxation(this)
            this.applyRelax = true;
        end

        %------------------------------------------------------------------
        function computePicardRelaxation(this, Xnew, XOld, DXOld)
            % Compute the current increment
            DX = Xnew - XOld;
            % Compute the generalized angle between successive increments
            gamma = acos((DX'*DXOld)/(norm(DX)*norm(DXOld)));
            % Update relaxation parameter
            if (gamma < pi/4.0)
                this.relax = min(this.relax * 2.0,1.0);
            elseif (gamma > pi/2.0)
                this.relax = this.relax / 2.0;
            end
        end


    end
    
end