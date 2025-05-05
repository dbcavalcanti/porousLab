%% EnrichedElement_H class
% This class extends the _RegularElement_H_ class to define a finite 
% element for single-phase fluid flow that incorporates enriched elements 
% to handle discontinuities. It provides methods to compute element data, 
% manage discontinuities, and calculate enriched degrees of freedom.
%
%% Methods
% * *elementData*: Computes the element data (stiffness matrix, damping 
%                  matrix, internal force vector, external force vector, 
%                  and derivative of internal force vector) based on 
%                  whether the element has discontinuities.
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
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef EnrichedElement_HM < RegularElement_HM   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_HM(node, elem, t, ...
                mat, intOrder, glu, glp, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress)
            this = this@RegularElement_HM(node, elem, t, ...
                mat, intOrder, glu, glp, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress);
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
           
           dfidu = zeros(this.ngle);
           if isempty(this.discontinuity)
               % Get the continuum contribution
               [Ke, Ce, fi, fe] = elementData@RegularElement_HM(this);
           else
               % Compute contribution of the continuum
               [Hcc, Hcj, Hjc, Hjj, Scc, Scj, Sjc, Sjj, fic, fij, fec, fej] = this.fillElementSubData();
               
               % Compute contribution of the discontinuity
               [Hdd, Sdd, Lcc, Lcj, Lcd, Ljc, Ljj, Ljd, Ldc, Ldj, Ldd, Tdc, Tdj] = getDiscontinuitiesData(this);
               
               % Condense the internal pressure dofs
               Hcc = Hcc + Tdc' * Hdd * Tdc;
               Hcj = Hcj + Tdc' * Hdd * Tdj;
               Hjc = Hjc + Tdj' * Hdd * Tdc;
               Hjj = Hjj + Tdj' * Hdd * Tdj;
               
               Scc = Scc + Tdc' * Sdd * Tdc;
               Scj = Scj + Tdc' * Sdd * Tdj;
               Sjc = Sjc + Tdj' * Sdd * Tdc;
               Sjj = Sjj + Tdj' * Sdd * Tdj;

               Lcc = Lcc + Tdc' * Ldd * Tdc - Lcd * Tdc - Tdc' * Ldc;
               Lcj = Lcj + Tdc' * Ldd * Tdj - Lcd * Tdj - Tdc' * Ldj;
               Ljc = Ljc + Tdj' * Ldd * Tdc - Ljd * Tdc - Tdj' * Ldc;
               Ljj = Ljj + Tdj' * Ldd * Tdj - Ljd * Tdj - Tdj' * Ldj;

               % Add contribution of the coupling matrices
               Hcc = Hcc + Lcc;
               Hcj = Hcj + Lcj;
               Hjc = Hjc + Ljc;
               Hjj = Hjj + Ljj;
               
               % Assemble the element matrices
               Ke = [ Hcc, Hcj;
                      Hjc, Hjj];
               
               Ce = [ Scc, Scj;
                      Sjc, Sjj];
               
               fi = [fic; fij];
               
               fe = [fec; fej];
           end
        end

        %------------------------------------------------------------------
        % Computes and assembles the sub-matrices and sub-vectors for an
        % enriched finite element
        function [Hcc, Hcj, Hjc, Hjj, Scc, Scj, Sjc, Sjj, fic, fij, fec, fej] = fillElementSubData(this)

            % Each discontinuity has a pressure jump dof
            nPJumpDofs = this.getNumberOfDiscontinuities();

            % Initialize fluid-flow sub-matrices
            Hcc = zeros(this.nglp , this.nglp);
            Hcj = zeros(this.nglp , nPJumpDofs);
            Hjc = zeros(nPJumpDofs, this.nglp);
            Hjj = zeros(nPJumpDofs, nPJumpDofs);

            % Initialize compressibity sub-matrices
            Scc = zeros(this.nglp , this.nglp);
            Scj = zeros(this.nglp , nPJumpDofs);
            Sjc = zeros(nPJumpDofs, this.nglp);
            Sjj = zeros(nPJumpDofs, nPJumpDofs);

            % Initialize external force vector
            fec = zeros(this.nglp, 1);
            fej = zeros(nPJumpDofs, 1);

            % Initialize internal force vector
            fic = zeros(this.nglp, 1);
            fij = zeros(nPJumpDofs, 1);
            
            % Vector of the nodal dofs
            pl = this.getNodalPressure();

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Cartesian coordinates of the integration point
                Xcar = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Enriched shape function matrix
                Npenr = this.enrichedShapeFncMtrx(Np, Xcar);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Compute the G matrix
                G = this.Gmatrix(Bp);

                % Pressure values at the integration point
                pIP = Np * pl;
        
                % Compute the permeability matrix
                kh = this.intPoint(i).constitutiveMdl.permeabilityTensor();

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibilityCoeff();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
        
                % Compute permeability sub-matrices
                Hcc = Hcc + Bp' * kh * Bp * c;
                Hcj = Hcj + Bp' * kh * G  * c;
                Hjc = Hjc + G'  * kh * Bp * c;
                Hjj = Hjj + G'  * kh * G  * c;

                % Compute compressibility matrices
                Scc = Scc + Np'    * comp * Np    * c;
                Scj = Scj + Np'    * comp * Npenr * c;
                Sjc = Sjc + Npenr' * comp * Np    * c;
                Sjj = Sjj + Npenr' * comp * Npenr * c;
                
                % Compute the gravity forces
                if (this.gravityOn)
                    fec = this.addGravityForces(fec,Bp,kh,pIP,c);
                end
            end

        end

        %------------------------------------------------------------------
        % Computes the stiffness matrix and force vector contributions
        % from discontinuities in the enriched element
        function [Hdd, Sdd, Lcc, Lcj, Lcd, Ljc, Ljj, Ljd, Ldc, Ldj, Ldd, Tdc, Tdj] = getDiscontinuitiesData(this)

            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nPJumpDofs        = nDiscontinuities;
            nPIntDofs         = nDiscontinuities * 2;

            % Initialize the output data 
            Hdd = zeros(nPIntDofs,nPIntDofs);
            Sdd = zeros(nPIntDofs,nPIntDofs);

            % Flow coupling matrices
            Lcc = zeros(this.nglp  , this.nglp);
            Lcj = zeros(this.nglp  , nPJumpDofs);
            Lcd = zeros(this.nglp  , nPIntDofs);
            Ljc = zeros(nPJumpDofs , this.nglp);
            Ljj = zeros(nPJumpDofs , nPJumpDofs);
            Ljd = zeros(nPJumpDofs , nPIntDofs);
            Ldc = zeros(nPIntDofs  , this.nglp);
            Ldj = zeros(nPIntDofs  , nPJumpDofs);
            Ldd = zeros(nPIntDofs  , nPIntDofs);

            % Condensation matrices
            Tdc = zeros(nPIntDofs,this.nglp);
            Tdj = zeros(nPIntDofs,nPJumpDofs);

            % Loop through the discontinuities
            for i = 1:nDiscontinuities

                % Internal pressure dofs associated with this discontinuity segment
                pint_dofs = 2*(i-1)+1 : 2*i;

                % Get the discontinuity data
                [Hddi, Sddi, Lcci, Lcji, Lcdi, Ljci, Ljji, Ljdi, Ldci, Ldji, Lddi,~,~,~] = this.discontinuity(i).elementData(this,i);

                % Assemble the contribution of this discontinuity
                Hdd(pint_dofs,pint_dofs) = Hdd(pint_dofs,pint_dofs) + Hddi;
                Sdd(pint_dofs,pint_dofs) = Sdd(pint_dofs,pint_dofs) + Sddi;
                Lcc = Lcc + Lcci;
                Lcj(:,i) = Lcj(:,i) + Lcji;
                Lcd(:,pint_dofs) = Lcd(:,pint_dofs) + Lcdi;
                Ljc(i,:) = Ljc(i,:) + Ljci;
                Ljj(i,i) = Ljj(i,i) + Ljji;
                Ljd(i,pint_dofs) = Ljd(i,pint_dofs) + Ljdi;
                Ldc(pint_dofs,:) = Ldc(pint_dofs,:) + Ldci;
                Ldj(pint_dofs,i) = Ldj(pint_dofs,i) + Ldji;
                Ldd(pint_dofs,pint_dofs) = Ldd(pint_dofs,pint_dofs) + Lddi;


                % Loop through the nodes of the discontinuity to fill the
                % condensation matrix
                for j = 1:2
                    X = this.discontinuity(i).node(j,:);
                    Xn = this.shape.coordCartesianToNatural(this.node,X);
                    Np = this.shape.shapeFncMtrx(Xn);
                    Tdc(pint_dofs(j),:) = Np;
                    Nenr_i = this.enrichedShapeFncValues(i, Np, X);
                    Tdj(pint_dofs(j),i) = 0.5 * (Nenr_i(2) + Nenr_i(3));
                end
                
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
        end

        %------------------------------------------------------------------
        % Adds a discontinuity segment to the element
        function addDiscontinuitySegment(this,dseg)
            this.discontinuity = [this.discontinuity; dseg];
        end

        %------------------------------------------------------------------
        % Adds the discontinuities dofs to the element dof vector
        function addEnrichmentToDofVector(this)
            nDiscontinuities = this.getNumberOfDiscontinuities();
            for i = 1:nDiscontinuities
                this.gle = [this.gle, this.discontinuity(i).dof];
            end
            this.ngle = length(this.gle);
        end

        %------------------------------------------------------------------
        % Compute the enrichment shape function matrix
        function Nenr = enrichedShapeFncMtrx(this, N, Xcar) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            Nenr = zeros(1,nDiscontinuities);  % Each discontinuity has a pressure jump dof
            for i = 1:nDiscontinuities  
                Nenr_i = this.enrichedShapeFncValues(i, N, Xcar);
                Nenr(i) = Nenr_i(1);
            end
        end

        %------------------------------------------------------------------
        % Compute the enrichment shape function values of a discontinuity
        function Nenr_i = enrichedShapeFncValues(this, id, N, Xcar) 
            phi = 0.0;
            h = this.discontinuity(id).heaviside(Xcar);
            for j = 1:this.nnd_el
                Xj = this.node(j,:);
                hj = this.discontinuity(id).heaviside(Xj);
                if (hj > 0.0)
                    phi = phi + N(:, j);
                end
            end
            Nenr_i = [h-phi;
                       -phi;        % Nbot
                      1.0-phi];     % Ntop
        end

        %------------------------------------------------------------------
        % Compute the gradient enrichment matrix
        function G = Gmatrix(this, Bu) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            G = zeros(2,nDiscontinuities);  % Each discontinuity has a pressure jump dof
            for i = 1:nDiscontinuities    
                for j = 1:this.nnd_el
                    Xj = this.node(j,:);
                    h = this.discontinuity(i).heaviside(Xj);
                    if (h > 0.0)
                        G(:,i) = G(:,i) - Bu(:, j);
                    end
                end
            end
        end
        %------------------------------------------------------------------
        
    end
end