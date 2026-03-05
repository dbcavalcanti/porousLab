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
classdef ECFractureElement_H < RegularElement_H   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ECFractureElement_H(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric)
            this = this@RegularElement_H(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric);
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
               [Ke, Ce, fi, fe] = elementData@RegularElement_H(this);
           else
               % Compute contribution of the continuum
               [Ke, Ce, fi, fe, dfidu] = this.fillECElementData();
           end
        end

        %------------------------------------------------------------------
        % Computes and assembles the sub-matrices and sub-vectors for an
        % enriched finite element
        function [Ke, Ce, fi, fe, dfidu] = fillECElementData(this)
      
            % Initialize the sub-matrices
            Ke    = zeros(this.nglp, this.nglp);
            Ce    = zeros(this.nglp, this.nglp);
            dfidu = zeros(this.nglp, this.nglp);

            % Initialize external force vector
            fe = zeros(this.nglp, 1);
            fi = zeros(this.nglp, 1);
            
            % Vector of the nodal dofs
            pl = this.getNodalPressure();

            % Initialize the volume of the element
            vol = 0.0;

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Pressure values at the integration point
                pIP = Np * pl;
        
                % Compute the permeability matrix
                kh = this.permeabilityTensor(i);

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibilityCoeff();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
        
                % Compute permeability sub-matrices
                Ke = Ke + Bp' * kh * Bp * c;

                % Compute compressibility matrices
                if ((this.massLumping) && (this.lumpStrategy == 1))
                    Ce = Ce + diag(comp*Np*c);
                elseif (this.massLumping == false)
                    Ce = Ce + Np' * comp * Np * c;
                end
                
                % Compute the gravity forces
                if (this.gravityOn)
                    fe = this.addGravityForces(fe,Bp,kh,pIP,c);
                end

                % Compute the element volume
                vol = vol + c;
            end

            % Compute the lumped mass matrix
            if ((this.massLumping) && (this.lumpStrategy == 2))
                Ce = lumpedCompressibilityMatrix(this, vol);
            end

        end

        %------------------------------------------------------------------
        function Kxy = permeabilityTensor(this, ip_id)
            Kxy = this.intPoint(ip_id).constitutiveMdl.intrinsicPermeabilityTensor();
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            % Add the longitudinal permeability
            for i = 1:nDiscontinuities
                % Get discontinuity properties
                kdl = this.discontinuity(i).intPoint(1).constitutiveMdl.longitudinalPermeability();
                ld  = this.discontinuity(i).ld();
                md  = this.discontinuity(i).tangentialVector();
                kt = this.discontinuity(i).mat.transversalPermeability;
                % Add it to the porous media tensor
                Kxy = Kxy + kt * kdl / ld * (md * md');
            end
            % Add the transversal permeability effect
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
            Kxy = this.intPoint(ip_id).constitutiveMdl.permeabilityTensor(Kxy);
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
            n = 0;  
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
            if isempty(dof_j)||isempty(dof_d)
                return
            end
            this.gle = [this.gle, dof_j, dof_d];
            this.ngle = length(this.gle);
        end
        %------------------------------------------------------------------
        
    end
end