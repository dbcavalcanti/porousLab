%% Anl_Nonlinear Class
%
% This class implements the solution of a nonlinear incremental-iterative 
% structural analysis.
%
% The code was adapted from the Anl_Nonlinear class from NUMA-TF 
% (https://gitlab.com/rafaelrangel/numa-tf, Accessed on February 1st, 2023)
% In the reference code to compute the stiffness matrix and the internal
% force vector are computed through different methods. Now both are
% computed using the same function. Another change that was done was
% related to the input variable to those methods, now the increment of the
% displacement vector associated to the iteration is used as an input.
% To consider the possibility of material nonlinearity, it was added a
% method to update the state variables after the convergence of the
% iterative process. It was also added the displacement control method.
%
classdef Anl_Nonlinear < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        method     = 0;   % flag for solution method
        adjustStep = 0;   % flag for type of increment size adjustment
        increment  = 0;   % initial increment of load ratio
        max_lratio = 0;   % limit value of load ratio
        max_step   = 0;   % maximum number of steps
        max_iter   = 0;   % maximum number of iterations in each step
        trg_iter   = 0;   % desired number of iterations in each step
        tol        = 0;   % numerical tolerance for convergence
        ctrlNode   = 0;   % control node (for displacement control method)
        ctrlDof    = 0;   % control dof (for displacement control method)
        incrSign   = 0;   % sign of the increment of displacement
        Uplot      = [];  % matrix of nodal displacement vectors of all steps/modes
        lbdplot    = [];  % vector of load ratios of all steps
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anl = Anl_Nonlinear(method,adjustStep,increment,...
                max_lratio,max_step,max_iter,trg_iter,tol)
            
            if nargin < 1, result = 0; end
            anl = anl@Anl('Nonlinear');

            % Default analysis configuration
            if nargin == 1
                anl.method     = 'LoadControl';
                anl.adjustStep = false;
                anl.increment  = 0.1;
                anl.max_lratio = 0.5;
                anl.max_step   = 40;
                anl.max_iter   = 10;
                anl.trg_iter   = 3;
                anl.tol        = 0.00001;
            
            else
            % User given analysis set up
                anl.method     = method;
                anl.adjustStep = adjustStep;
                anl.increment  = increment;
                anl.max_lratio = max_lratio;
                anl.max_step   = max_step;
                anl.max_iter   = max_iter;
                anl.trg_iter   = trg_iter;
                anl.tol        = tol;
            end

            % For the displacement control method
            anl.ctrlDof = result.dof;
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function status = process(anl,mdl)

            % Initialize the model object
            mdl.preComputations();

            % Initialize results
            anl.lbdplot = zeros(anl.max_step+1,1);
            anl.Uplot   = zeros(anl.max_step+1,1);
            
            % Initialize data for first step
            step  = 0;  % step number
            lbd   = 0;  % total load ratio (lambda)
            sign  = 1;  % sign of predicted increment of load ratio
            
            % Initialize vector of total nodal displacements
            U    = mdl.U;

            % Initialize vector of total increment displacement
            D_U = zeros(mdl.ndof,1);
            
            %==========================================================================
            % Start incremental process
            while (step < anl.max_step)
                step = step + 1;
                
                % Tangent stiffness matrix
                [K, ~, ~, Fref] = mdl.globalMatrices(U);

                % Tangent increment of displacements for predicted solution
                d_Up0 = anl.solveSystem(mdl,K,Fref,U);
                
                if (step == 1)
                    % Initial increment of load ratio for predicted solution
                    if strcmp(anl.method,'DisplacementControl')
                        d_lbd0 = anl.predictedIncrement(anl,mdl,sign,1,1,0.0,0.0,D_U,d_Up0,Fref);
                    else
                        d_lbd0 = anl.increment;
                    end
                    
                    % Set previous tangent increment of displacements as current increment
                    d_Up0_old = d_Up0;
                    
                    % Store squared value of the norm of tangent increment of displacements
                    n2 = d_Up0(mdl.doffree)'*d_Up0(mdl.doffree);
                else
                    % Generalized Stiffness Parameter
                    GSP = n2/(d_Up0(mdl.doffree)'*d_Up0_old(mdl.doffree));
                    
                    % Adjust increment sign
                    if (GSP < 0)
                        sign = -sign;
                    end
                    
                    % Adjustment factor of increment size
                    if (anl.adjustStep == false)
                        J = 1;
                    elseif (anl.adjustStep == true)
                        J = sqrt(anl.trg_iter/iter);
                    end
                    
                    % Predicted increment of load ratio
                    d_lbd0 = anl.predictedIncrement(anl,mdl,sign,J,GSP,D_lbd,d_lbd0,D_U,d_Up0,Fref);
                end
                
                % Limit increment of load ratio to make total load ratio smaller than maximum value
                if ((anl.max_lratio > 0.0 && lbd + d_lbd0 > anl.max_lratio) ||...
                    (anl.max_lratio < 0.0 && lbd + d_lbd0 < anl.max_lratio))
                    d_lbd0 = anl.max_lratio - lbd;
                end
                
                % Increments of load ratio and displacements for predicted solution
                d_lbd = d_lbd0;
                d_U0  = d_lbd0 * d_Up0;
                d_U   = d_U0;
                
                % Initialize incremental values of load ratio and displacements for current step
                D_lbd = d_lbd;
                D_U   = d_U;
                
                % Update total values of load ratio and displacements
                lbd = lbd + d_lbd;
                U   = U   + d_U;
                
                %----------------------------------------------------------------------
                % Start iterative process
                iter = 1;
                conv = 0;
                while (conv == 0 && iter <= anl.max_iter)
                    
                    % Vector of external and internal forces
                    Fext = lbd * Fref;
                    [K, ~, Fint] = mdl.globalMatrices(U);
                    
                    % Vector of unbalanced forces
                    R = Fext - Fint;
                    
                    % Check convergence
                    unbNorm = norm(R(mdl.doffree));
                    forNorm = norm(Fref(mdl.doffree));
                    conv = (unbNorm == 0 || forNorm == 0 || unbNorm/forNorm < anl.tol);
                    if conv == 1
                        break;
                    end
                    
                    % Start/keep corrector phase
                    iter = iter + 1;
                       
                    % Tangent and residual increments of displacements
                    d_Up = anl.solveSystem(mdl,K,Fref);
                    d_Ur = anl.solveSystem(mdl,K,R);
                    
                    % Corrected increment of load ratio
                    d_lbd = anl.correctedIncrement(anl,mdl,d_lbd0,D_lbd,d_Up0_old,d_U0,d_Up,d_Ur,D_U,Fref,R);
                    if (~isreal(d_lbd))
                        conv = -1;
                        break;
                    end
                    
                    % Corrected increment of displacements
                    d_U = d_lbd * d_Up + d_Ur;
                    
                    % Increments of load ratio and displacements for current step
                    D_lbd = D_lbd + d_lbd;
                    D_U   = D_U   + d_U;
                    
                    % Total values of load ratio and displacements
                    lbd = lbd + d_lbd;
                    U   = U   + d_U;
                end
                %----------------------------------------------------------------------
                % Check for convergence fail or complex value of increment
                if (conv == 0)
                    status = (step > 1);
                    fprintf('Status: Convergence not achieved!\n');
                    anl.plotCurves();
                    return;
                elseif (conv == -1)
                    status = (step > 1);
                    fprintf('Status: Unable to compute load increment!\n');
                    return;
                end

                fprintf('Step:%d | Iter:%d | ratio:%.2f\n',step,iter,lbd);

                % Update the state variables
                mdl.updateStateVar();
                
                % Store step results
                anl.lbdplot(step+1) = lbd;
                anl.Uplot(step+1) = U(res.dof);
                
                % Store predicted tangent increment of displacements for next step
                if (step ~= 1)
                    d_Up0_old = d_Up0;
                end
                
                % Check if maximum load ratio was reached
                if ((anl.max_lratio >= 0 && lbd >= 0.999*anl.max_lratio) ||...
                    (anl.max_lratio <= 0 && lbd <= 0.999*anl.max_lratio))
                    break;
                end
            end
            %==========================================================================
            
            % Clean unused steps
            if (step < anl.max_step)
                anl.lbdplot = anl.lbdplot(1:step+1);
                anl.Uplot = anl.Uplot(1:step+1);
            end
            anl.plotCurves();

            mdl.U = U;

        end

        %------------------------------------------------------------------
        function plotCurves(this)
            figure
            hold on
            box on, grid on, axis on
            plot(this.Uplot, this.lbdplot, 'o-k');
            xlabel('Displacement (mm)','Interpreter','latex')
            ylabel('Load factor','Interpreter','latex')
            xaxisproperties= get(gca, 'XAxis');
            xaxisproperties.TickLabelInterpreter = 'latex'; 
            yaxisproperties= get(gca, 'YAxis');
            yaxisproperties.TickLabelInterpreter = 'latex';   
            set(gca,'FontSize',14);
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Compute inrement of load ratio for the predicted solution
        % (first iteration)
        function d_lbd0 = predictedIncrement(anl,mdl,sign,J,GSP,D_lbd,d_lbd0,D_U,d_Up0,Pref)
            % Extract free d.o.f. components
            Pref  = Pref(mdl.doffree);
            D_U   = D_U(mdl.doffree);
            d_Up0 = d_Up0(mdl.doffree);
            
                
            % LCM: Load Increment
            if strcmp(anl.method,'LoadControl')
                d_lbd0 = J * abs(d_lbd0);

            % DCM: Displacement Increment
            elseif strcmp(anl.method,'DisplacementControl')
                d_lbd0 = J * sign * anl.increment / d_Up0(anl.ctrlDof);
                return  % The change of the sign must not be applied
                
            % WCM: Work Increment
            elseif strcmp(anl.method,'WorkControl')
                d_lbd0 = J * sqrt(abs((D_lbd*Pref'*D_U)/(Pref'*d_Up0)));
                
            % ALCM_FNP: Cylindrical Arc-Length Increment
            elseif strcmp(anl.method,'ArcLengthFNPControl')
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % ALCM_UNP: Cylindrical Arc-Length Increment
            elseif strcmp(anl.method,'ArcLengthUNPControl')
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % ALCM_CYL: Cylindrical Arc-Length Increment
            elseif strcmp(anl.method,'ArcLengthCylControl')
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % ALCM_SPH: Spherical Arc-Length Increment
            elseif strcmp(anl.method,'ArcLengthSPHControl')
                d_lbd0 = J * sqrt((D_U'*D_U + D_lbd^2*(Pref'*Pref)) / (d_Up0'*d_Up0 + Pref'*Pref));
                
            % MNCM: Cylindrical Arc-Length Increment
            elseif strcmp(anl.method,'MinimumNorm')
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % ORCM: Cylindrical Arc-Length Increment
            elseif strcmp(anl.method,'OrthogonalResidual')
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % GDCM: GSP criteria
            elseif strcmp(anl.method,'GeneralizedDisplacement')
                d_lbd0 = J * sqrt(abs(GSP)) * anl.increment;
            end
            
            % Apply increment sign
            d_lbd0 = sign * d_lbd0;
        end
        
        %--------------------------------------------------------------------------
        % Compute inrement of load ratio for the corrected solutions
        % (iterations to correct predicted solution).
        function d_lbd = correctedIncrement(anl,mdl,d_lbd0,D_lbd,d_Up0,d_U0,d_Up,d_Ur,D_U,Pref,R)
            % Extract free d.o.f. components
            d_Up0 = d_Up0(mdl.doffree);
            d_U0  = d_U0(mdl.doffree);
            d_Up  = d_Up(mdl.doffree);
            d_Ur  = d_Ur(mdl.doffree);
            D_U   = D_U(mdl.doffree);
            Pref  = Pref(mdl.doffree);
            R     = R(mdl.doffree);
            
            % LCM
            if strcmp(anl.method,'LoadControl')
                d_lbd = 0;
            
            % DCM
            elseif strcmp(anl.method,'DisplacementControl')
                d_lbd = -d_Ur(anl.ctrlDof)/d_Up(anl.ctrlDof);
                
            % WCM
            elseif strcmp(anl.method,'WorkControl')
                d_lbd = -(Pref'*d_Ur)/(Pref'*d_Up);
                
            % ALCM_FNP
            elseif strcmp(anl.method,'ArcLengthFNPControl')
                d_lbd = -(d_Ur'*d_U0)/(d_Up'*d_U0 + d_lbd0*(Pref'*Pref));
                
            % ALCM_UNP
            elseif strcmp(anl.method,'ArcLengthUNPControl')
                d_lbd = -(d_Ur'*D_U)/(d_Up'*D_U + D_lbd*(Pref'*Pref));
                
            % ALCM_CYL
            elseif strcmp(anl.method,'ArcLengthCylControl')
                a = d_Up'*d_Up;
                b = d_Up'*(d_Ur + D_U);
                c = d_Ur'*(d_Ur + 2*D_U);
                s = sign(D_U'*d_Up);
                
                d_lbd = -b/a + s*sqrt((b/a)^2 - c/a);
                
            % ALCM_SPH
            elseif strcmp(anl.method,'ArcLengthSPHControl')
                a = d_Up'*d_Up + Pref'*Pref;
                b = d_Up'*(d_Ur + D_U) + D_lbd*(Pref'*Pref);
                c = d_Ur'*(d_Ur + 2*D_U);
                s = sign(D_U'*d_Up);
                
                d_lbd = -b/a + s*sqrt((b/a)^2 - c/a);
                
            % MNCM
            elseif strcmp(anl.method,'MinimumNorm')
                d_lbd = -(d_Up'*d_Ur)/(d_Up'*d_Up);
                
            % ORCM
            elseif strcmp(anl.method,'OrthogonalResidual')
                d_lbd = -(R'*D_U)/(Pref'*D_U);
                
            % GDCM
            elseif strcmp(anl.method,'GeneralizedDisplacement')
                d_lbd = -(d_Up0'*d_Ur)/(d_Up0'*d_Up);
            end
        end
        
    end
end