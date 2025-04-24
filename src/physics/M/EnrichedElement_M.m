%% EnrichedElement_M Class
% This class defines an enriched mechanical finite element that extends 
% the functionality of the _RegularElement_M_ class. It incorporates 
% additional degrees of freedom to handle discontinuities and enrichment 
% modes such as stretching and relative rotation. 
%
%% Methods
% * *elementData*: Computes the element data (stiffness matrix, damping 
%                  matrix, internal force vector, external force vector, 
%                  and derivative of internal force vector) based on 
%                  whether the element has discontinuities.
% * *enrichedElementData*: Computes the enriched element data using static 
%                          condensation of the enrichment degrees of 
%                          freedom.
% * *solveLocalEq*: Solves the local equilibrium equation iteratively 
%                   using the Newton-Raphson method.
% * *fillElementSubData*: Fills the sub-matrices and vectors for the 
%                         enriched element based on the current enriched 
%                         degrees of freedom.
% * *getDiscontinuitiesData*: Computes the stiffness matrix and force 
%                             vector contributions from the 
%                             discontinuities.
% * *getNumberEnrichedDofs*: Returns the total number of enriched degrees 
%                            of freedom.
% * *getNumberOfDiscontinuities*: Returns the number of discontinuities 
%                                 associated with the element.
% * *getNumberOfDofPerDiscontinuity*: Returns the number of degrees of 
%                                     freedom per discontinuity, 
%                                     considering the enabled enrichment 
%                                     modes.
% * *addDiscontinuitySegment*: Adds a discontinuity segment to the element.
% * *kinematicEnrichment*: Computes the kinematic enrichment matrix for 
%                          the element based on the discontinuities and 
%                          enrichment modes.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00: Initial version (January 2024).
%
%% Class definition
classdef EnrichedElement_M < RegularElement_M    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
        addStretchingMode = false;
        addRelRotationMode = false;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_M(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress, ...
                addRelRotationMode,addStretchingMode)
            this = this@RegularElement_M(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress);
            this.addStretchingMode  = addStretchingMode;
            this.addRelRotationMode = addRelRotationMode;
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Computes the element data for the current element based on wether
        % the element contains a discontinuity or not.
        % 
        % Outputs:
        %   Ke    - Element stiffness matrix.
        %   Ce    - Element damping matrix.
        %   fi    - Internal force vector.
        %   fe    - External force vector.
        %   dfidu - Derivative of internal force with respect to 
        %           displacement.
        function [Ke, Ce, fi, fe, dfidu] = elementData(this)

           if isempty(this.discontinuity)
               [Ke, Ce, fi, fe, dfidu] = elementData@RegularElement_M(this);
           else
               [Ke, Ce, fi, fe, dfidu] = enrichedElementData(this);
           end
            
        end

        %------------------------------------------------------------------
        % Computes the enriched element data for the finite element.
        % 
        % Outputs:
        %   Ke    - Element stiffness matrix.
        %   Ce    - Element damping matrix.
        %   fi    - Internal force vector.
        %   fe    - External force vector.
        %   dfidu - Derivative of internal force with respect to 
        %           displacement.
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
        %  Solves the local equilibrium equation for an enriched element
        function [Kuu, Kua, Kau, Kaa, fiu, fia, fe] = solveLocalEq(this)

            % Initialize the enrichment dofs vector
            nEnrDofs = this.getNumberEnrichedDofs();
            ae = zeros(nEnrDofs,1);

            % Define Newton-Raphson parameter
            conv    = false;
            tol     = 1.0e-5;
            maxIter = 10;

            % Iterative solution of the local equilibrium equation
            for i = 1:maxIter

                % Compute sub-matrices
                [Kuu, Kua, Kau, Kaa, fiu, fia, fe] = this.fillElementSubData(ae);

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
        % Computes and assembles the sub-matrices and sub-vectors for an
        % enriched finite element
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

            % Initialize external force vector
            fe = zeros(this.nglu, 1);

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

                % Get the static enriched matrix (TEMP)
                Gv = this.kinematicEnrichment(Bu);

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
                Kuu = Kuu + Bu' * Duu * Bu * c;
                Kua = Kua + Bu' * Duu * Gr * c;
                Kau = Kau + Gv' * Duu * Bu * c;
                Kaa = Kaa + Gv' * Duu * Gr * c;

                % Internal force vector
                fiu = fiu + Bu' * stress * c;
                fia = fia + Gv' * stress * c;
                
                % Compute the gravity forces
                if (this.gravityOn)
                    fe = this.addGravityForces(fe,this.intPoint(i).X,c);
                end
            end

            % Add the contribution from the discontinuities
            [Kd,fd] = this.getDiscontinuitiesData(ae);
            fia = fia + fd;
            Kaa = Kaa + Kd;
        end

        %------------------------------------------------------------------
        % Computes the stiffness matrix and force vector contributions
        % from discontinuities in the enriched element
        function [Kd,fd] = getDiscontinuitiesData(this,ae)

            nEnrDofs          = this.getNumberEnrichedDofs();
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();

            % Initialize the output data 
            Kd = zeros(nEnrDofs,nEnrDofs);
            fd = zeros(nEnrDofs,1);

            % Loop through the discontinuities
            for i = 1:nDiscontinuities

                % Dofs associated with this discontinuity segment
                dofs = nDofDiscontinuity*(i-1)+1 : nDofDiscontinuity*i;

                % Get the discontinuity data
                [Kdi,~,fdi,~,~] = this.discontinuity(i).elementData(ae(dofs));

                % Assemble the contribution of this discontinuity
                Kd(dofs,dofs) = Kdi;
                fd(dofs)      = fdi;
                
            end
        end

        %------------------------------------------------------------------
        % Gets the number of enriched degrees of freedom
        function nEnrDof = getNumberEnrichedDofs(this)
            nEnrDof = this.getNumberOfDiscontinuities();
            nEnrDof = nEnrDof * this.getNumberOfDofPerDiscontinuity();
        end

        %------------------------------------------------------------------
        % Obtain the number of discontinuities
        function n = getNumberOfDiscontinuities(this)
            n = size(this.discontinuity,1);
        end

        %------------------------------------------------------------------
        % Calculates the number of degrees of freedom per discontinuity for
        % the enriched element
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
        % Adds a discontinuity segment to the element
        function addDiscontinuitySegment(this,dseg)
            this.discontinuity = [this.discontinuity; dseg];
        end

        %------------------------------------------------------------------
        % Computed the kinematic enrichment matrix for an enriched finite
        % element.
        % It calculated the enrichment matrix by considering the
        % contributions of discontinuities in the element. The enrichment
        % includes translation, stretching and relative rotation modes
        % depending on the configuration
        function Gc = kinematicEnrichment(this, Bu) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();
            Gc = zeros(4,nDofDiscontinuity * nDiscontinuities);
            for i = 1:nDiscontinuities    
                Gci = zeros(4,nDofDiscontinuity);
                % Get the discontinuity orientation vectors
                m = this.discontinuity(i).tangentialVector();
                n = this.discontinuity(i).normalVector();
                % Get the discontinuity reference point
                Xr = this.discontinuity(i).referencePoint();
                for j = 1:this.nnd_el
                    Xj = this.node(j,:);
                    h = this.discontinuity(i).heaviside(Xj);
                    if (h > 0.0)
                        % Columns of the B-matrix associated with this node
                        Buj = Bu(:, 2*(j-1) + 1 : 2*j);
                        % Add translation modes
                        Gci(:,1) = Gci(:,1) - Buj * m;
                        Gci(:,2) = Gci(:,2) - Buj * n;
                        % Add stretching mode
                        c = 3;
                        if this.addStretchingMode
                            Gci(:,c) = Gci(:,c) - Buj * (m * m') * (Xj' - Xr');
                            c = c + 1;
                        end
                        % Add relative rotation mode
                        if this.addRelRotationMode
                            mn = (n * m') - (m * n');
                            Gci(:,c) = Gci(:,c) - Buj * mn * (Xj' - Xr');
                        end
                    end
                end
                % Assemble the matrix associated with discontinuity i
                cols = nDofDiscontinuity*(i-1)+1 : nDofDiscontinuity*i;
                Gc(:,cols) = Gci;
            end
        end
        %------------------------------------------------------------------
    end
end