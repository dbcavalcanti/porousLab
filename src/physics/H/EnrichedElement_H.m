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
% * *pressureField*: Computes the matrix pressure, including the pressure
%                    jumps across embedded discontinuities.
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef EnrichedElement_H < RegularElement_H   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
        condenseInternalPressure = true;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_H(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,condenseDiscontinuityPressure)
            this = this@RegularElement_H(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric);
            this.condenseInternalPressure = condenseDiscontinuityPressure;
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Computes the element data for the current element based on whether
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
               [Ke, Ce, fi, fe] = elementData@RegularElement_H(this);
           else
               % Compute contribution of the continuum
               [Hcc, Hcj, Hjc, Hjj, Scc, Scj, Sjc, Sjj, fic, fij, fec, fej] = this.fillElementSubData();
               
               % Compute contribution of the discontinuity
               [Hdd, Sdd, Lcc, Lcj, Lcd, Ljc, Ljj, Ljd, Ldc, Ldj, Ldd, Tdc, Tdj] = getDiscontinuitiesData(this);

               % Zero auxiliary matrices
               Ocd = zeros(this.nglp, size(Hdd,1));
               Ojd = zeros(size(Hjj,1), size(Hdd,1));
               od  = zeros(size(Hdd,1),1);

               % Assemble the flow matrix
               if this.condenseInternalPressure
                   D = [Tdc, Tdj];
                   Ke = [Hcc, Hcj; Hjc, Hjj] + D' * Hdd * D;

                   % With endpoint quadrature and equal leakoff on both
                   % faces, the constrained internal pressure leaves only
                   % the jump penalty, integral(leakoff/2 ds). Assemble it
                   % directly to avoid cancellation of large coupling terms.
                   for i = 1:this.getNumberOfDiscontinuities()
                       pint_dofs = 2*(i-1)+1 : 2*i;
                       jump_dof = this.nglp + i;
                       penalty = sum(sum(Ldd(pint_dofs,pint_dofs))) / 4;
                       Ke(jump_dof,jump_dof) = Ke(jump_dof,jump_dof) + penalty;
                   end
               else
                   Ke = [ Hcc+Lcc, Hcj+Lcj,   -Lcd;
                          Hjc+Ljc, Hjj+Ljj,   -Ljd;
                             -Ldc,    -Ldj, Hdd+Ldd];
               end

               Ce = [ Scc,  Scj,  Ocd;
                      Sjc,  Sjj,  Ojd;
                      Ocd', Ojd', Sdd];

               fi = [fic; fij; od];
               
               fe = [fec; fej; od];
               
               % Condense the internal pressure dofs
               if this.condenseInternalPressure
                   T = this.condensationMatrix(Tdc, Tdj);
                   Ce = T' * Ce * T;
                   fi = T' * fi;
                   fe = T' * fe;
               end
               
               
           end
        end

        %------------------------------------------------------------------
        % Computes the matrix pressure, including all discontinuity jumps.
        % X is given in global cartesian coordinates; ue is optional.
        % Internal fracture pressures are separate from this matrix field.
        function p = pressureField(this,X,ue)
            if nargin > 2, this.ue = ue; end

            % Natural coordinates and regular shape functions
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            Np = this.shape.shapeFncMtrx(Xn);

            % Nodal contribution
            p = Np * this.getNodalPressure();

            % Jump dofs follow the nodal pressures in both the condensed
            % and explicit internal-pressure formulations.
            nDiscontinuities = this.getNumberOfDiscontinuities();
            if nDiscontinuities > 0
                Nenr = this.enrichedShapeFncMtrx(Np,X);
                pj = this.ue(this.nglp+(1:nDiscontinuities));
                p = p + Nenr * pj;
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
                kh = permeabilityTensor(this, i);

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

        function Kxy = permeabilityTensor(this, ip_id)
            Kxy = this.intPoint(ip_id).constitutiveMdl.permeabilityTensor();
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            % Loop through the discontinuities
            for i = 1:nDiscontinuities
                R = this.discontinuity(i).rotationFromGlobalToLocal();
                % Rotate the permeability tensor
                Kmn = R * Kxy * R';
                % Apply transversal permeability coefficient
                kt = this.discontinuity(i).mat.transversalPermeability;
                kt = min(max(kt, 0.0), 1.0);
                Kmn(2,2) = kt * Kmn(2,2);
                % Rotate back to the cartesian system
                Kxy = R' * Kmn * R;
            end
        end

        %------------------------------------------------------------------
        % Internal pressure dof condensation matrix
        function T = condensationMatrix(this, Tdc, Tdj)
            % Number of dofs of each type
            ndof_c = this.nglp;
            ndof_j = size(Tdj, 2);
            % Auxiliary matrices
            Icc = eye(ndof_c);
            Ijj = eye(ndof_j);
            Ocj = zeros(ndof_c, ndof_j);
            % Condensation matrix
            T = [ Icc , Ocj;
                  Ocj', Ijj;
                  Tdc, Tdj];
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
            n = 1;
            if ~this.condenseInternalPressure
                n = n + 2;
            end
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
            dof_j = []; dof_d = [];
            for i = 1:nDiscontinuities
                dof_j = [dof_j, this.discontinuity(i).dof(1)];
                dof_d = [dof_d, this.discontinuity(i).dof(2:end)];
            end
            if isempty(dof_j)
                return
            end
            if this.condenseInternalPressure
                this.gle = [this.glp, dof_j];
            else
                this.gle = [this.glp, dof_j, dof_d];
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
