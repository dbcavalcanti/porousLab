%% EnrichedElement_M class
%
% This class defines a mechanical finite element 
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef EnrichedElement_M < RegularElement_M    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
        nDiscontinuities = 0;
        addStretchingMode = false;
        addRelRotationMode = false;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_M(type, node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress)
            this = this@RegularElement_M(type, node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % This function assembles the element matrices and vectors 
        %
        % Output:
        %   Ke : element "stiffness" matrix
        %   Ce : element "damping" matrix
        %   fe : element "internal force" vector
        %
        function [Ke, Ce, fi, fe, dfidu] = elementData(this)

           if isempty(this.discontinuity)
               [Ke, Ce, fi, fe, dfidu] = elementData@RegularElement_M(this);
           else
               [Ke, Ce, fi, fe, dfidu] = enrichedElementData(this);
           end
            
        end

        %------------------------------------------------------------------
        function [Ke, Ce, fi, fe, dfidu] = enrichedElementData(this)

            % Initialize the matrices and vectors that will be returned
            Ce    = zeros(this.nglu, this.nglu);
            dfidu = zeros(this.nglu, this.nglu);

            % Compute the sub-matrices
            [Kuu, Kua, Kau, Kaa, fiu, fia, fe] = this.solveLocalEq();

            % Static condensation of the enrichment dofs
            Ke = Kuu - Kua * (Kaa\Kau);
            fi = fiu - Kua * (Kaa\fia);

        end

        %------------------------------------------------------------------
        function [Kuu, Kua, Kau, Kaa, fiu, fia, fe] = solveLocalEq(this)

            % Initialize the enrichment dofs vector
            nEnrDofs = this.getNumberEnrichedDofs();
            ae = zeros(nEnrDofs,1);

            % Define Newton-Raphson parameter
            conv = false;
            tol = 1.0e-5;

            % Iterative solution of the local equilibrium equation
            for i = 1:10

                % Compute sub-matrices
                [Kuu, Kua, Kau, Kaa, fiu, fia, fe] = this.fillElementSubData(this);

                % Check convergence
                if (norm(fia) < tol)
                    conv = true;
                    break
                end

                % Solve iterative equation
                dae = -Kaa\fia;

                % Update enriched dofs
                ae = ae + dae;
                    
            end

            % Stop analysis process if solution did not converge
            if (conv == false)
                disp('LOCAL EQUILIBRIUM EQUATION FAIL TO CONVERGE');
                error('Local equilibrium equation did not converge');
            end

        end

        %------------------------------------------------------------------
        function [Kuu, Kua, Kau, Kaa, fiu, fia, fe] = fillElementSubData(this,ae)
            % Get the number of dofs of each type
            nEnrDofs = this.getNumberEnrichedDofs();
            nRegDofs = this.nglu;

            % Initialize the sub-matrices
            Kuu = zeros(nRegDofs,nRegDofs);
            Kua = zeros(nRegDofs,nEnrDofs);
            Kau = zeros(nEnrDofs,nRegDofs);
            Kaa = zeros(nEnrDofs,nEnrDofs);

            % Initialize the sub-vectors
            fiu = zeros(nRegDofs,1);
            fia = zeros(nEnrDofs,1);

            % Vector of the nodal dofs
            u  = this.getNodalDisplacement();

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints
               
                % Compute the B matrix at the int. point and the detJ
                [dNdx, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.shape.BMatrix(dNdx);

                % Get kinematic enriched matrix
                Gr = this.kinematicEnrichment(Bu);

                % Compute the strain vector
                this.intPoint(i).strain = Bu * u + Gr * ae;

                % Compute the stress vector and the constitutive matrix
                [stress,Duu] = this.intPoint(i).mechanicalLaw();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
                
                % Compute the stiffness sub-matrix
                Ke = Ke + Bu' * Duu * Bu * c;

                % Internal force vector
                fi = fi + Bu' * stress * c;
                
                % Compute the gravity forces
                if (this.mat.porousMedia.gravityOn)
                    fe = this.addGravityForces(fe,this.intPoint(i).X,c);
                end
            end 
        end

        %------------------------------------------------------------------
        function nEnrDof = getNumberEnrichedDofs(this)
            nEnrDof = this.getNumberOfDiscontinuities();
            nEnrDof = nEnrDof * this.getNumberOfDofPerDiscontinuity();
        end

        %------------------------------------------------------------------
        function n = getNumberOfDiscontinuities(this)
            n = size(this.discontinuity,1);
        end

        %------------------------------------------------------------------
        function n = getNumberOfDofPerDiscontinuity(this)
            n = 2;  
            if this.addStretchingMode
                n = n + 1;
            end
            if this.addRelRotationMode
                n = n + 1;
            end
        end

        %------------------------------------------------------------------
        function addDiscontinuitySegment(this,dseg)
            this.discontinuity = [this.discontinuity; dseg];
        end

        %------------------------------------------------------------------
        function kinematicEnrichment(this, Bu)

        end


    end
end