%% EnrichedElement class
%
% This class defines an enriched finite element
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: December 2022
%%%
%
%% Class definition
classdef EnrichedElementSmoothed < EnrichedElement
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        neighbors = [];
        nglu_neighbor = [];
        gle_smoothed = [];
        displJump = [];
        Kwd = [];
    end
    %% Constructor method
    methods
        function this = EnrichedElementSmoothed(type, node, elem, anm, t, matModel, mat, intOrder, glu, glp, fracture, subDivInt, stretch, enrVar, jumpOrder, staticCondensation, fractIntersectionFlag, nFractIntersectionDof, fractIntersectionId, EAS)
            this = this@EnrichedElement(type, node, elem, anm, t, matModel, mat, intOrder, glu, glp, fracture, subDivInt, stretch, enrVar, jumpOrder, staticCondensation, fractIntersectionFlag, nFractIntersectionDof, fractIntersectionId, EAS);
        end
    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        function fillDofData(this)

            % Get the regular dofs of the neighbor elements
            gluN = []; this.nglu_neighbor = zeros(length(this.neighbors),1);
            for i = 1:length(this.neighbors)
                if this.staticCondensation.DisplJump
                    gluN = [gluN, this.neighbors(i).glu];
                    this.nglu_neighbor(i) = length(this.neighbors(i).glu);
                else
                    gluN = [gluN, this.neighbors(i).gluenr];
                    this.nglu_neighbor(i) = this.neighbors(i).ngluenr;
                end
            end

            % Assemble the element dofs
            this.gle          = [this.glu, this.gluenr, this.glp , this.glpenrDp , this.glpenrPf];
            this.gle_smoothed = [this.glu, this.gluenr, gluN , this.glp , this.glpenrDp , this.glpenrPf];
            this.ngle         = length(this.gle);
        end

        %------------------------------------------------------------------
        % This function assembles the element matrices and vectors 
        % Input:
        %   dUe: vector with increment of the nodal dof vector
        %        associated with the element.
        %
        % Output:
        %   Ke : element "stiffness" matrix
        %   Ce : element "damping" matrix
        %   fe : element "internal force" vector
        %
        %   Ke = [ K     -Q  ]
        %        [ 0     H+L ]
        %
        %   Ce = [  0    0 ]
        %        [ Q^T   S ]
        %
        function [Ke, Ce, fe] = elementData(this)

            % Compute the smoothed stiffness matrix
            [K, fi] = this.elementSmoothedStiffnessMatrix();

            % Number of dofs associated with the neighbors
            ngluNeighbors = sum(this.nglu_neighbor(:));
            
            % Initialize auxiliary zero matrices
            Ouu = zeros(this.nglutot + ngluNeighbors, this.nglutot + ngluNeighbors);
            Opu = zeros(this.nglptot, this.nglutot + ngluNeighbors);
            Opp = zeros(this.nglptot, this.nglptot);

            % Local degrees of freedom
            uDofs = 1:(this.nglutot + ngluNeighbors);

            % Initialize element internal force vector
            fe = zeros(length(this.gle_smoothed), 1);

            % Assemble internal force vector
            fe(uDofs) = fi(uDofs);

            % Assemble the matrices
            Ke = [  K  , Opu';
                   Opu , Opp ];

            Ce = [ Ouu ,  Opu';
                   Opu ,  Opp ];
    
            
        end

        %------------------------------------------------------------------
        % Compute the smoothed stiffness matrix
        function [K,fi] = elementSmoothedStiffnessMatrix(this)

            nNeighbors = length(this.neighbors);
            nglw_neighbors = sum(this.nglu_neighbor(:)); 

            % Initialize
            Kc = zeros(this.nglutot+nglw_neighbors, this.nglutot+nglw_neighbors);  % Stiffness matrix 
            fi = zeros(this.nglutot+nglw_neighbors, 1);

            % Get element displacement jumps
            UeEnr = this.getDisplacementJumps();
            Ai = this.getDomainArea();
            At = Ai;
            Dist = zeros(1,nNeighbors);
            for j = 1:nNeighbors
                if nNeighbors > 1
                    % factor = 1.0/length(this.neighbors(j).neighbors);
                        % factor = 1.0/nNeighbors;
                        factor = 0.5;
                else
                    factor = 1.0;
                end
                At    = At + factor*this.neighbors(j).getDomainArea();
                UeEnr = [UeEnr; this.neighbors(j).getDisplacementJumps()];
            end

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints

                % Compute the B matrix at the int. point and the detJ
                [dNdx, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.shape.BMatrix(dNdx);

                % Compute the enriched B-matrix
                [~, BuEnr] = this.enhancedCompatibilityMtrx(dNdx,this.intPoint(i).X);

                Gkson = this.enhancedStaticMatrix(this.fracture{1}(1));

                % Compute the smoothed enhanced strain 
                BuEnr = BuEnr * Ai;
                Gkson = Gkson * Ai;
                for j = 1:nNeighbors
                    % Get the neighbor element data
                    An = this.neighbors(j).getDomainArea();
                    if nNeighbors > 1
                        % factor = 1.0/length(this.neighbors(j).neighbors);
                        % factor = 1.0/nNeighbors;
                        factor = 0.5;
                    else
                        factor = 1.0;
                    end
                    BuEnr = [BuEnr, factor*An*this.neighbors(j).getBEnrMatrix()];
                    Gkson = [Gkson, factor*An*this.neighbors(j).enhancedStaticMatrix(this.neighbors(j).fracture{1}(1));];
                end
                BuEnr = BuEnr / At;
                Gkson = Gkson / At;
                
                % Compute the increment of the strain vector
                strain = Bu*this.ue(1:this.nglu) + BuEnr * UeEnr;

                % Compute the stress vector and the constitutive matrix
                [stress,Du] = this.intPoint(i).constitutiveModel(strain);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;

                % Construct augmented matrices
                Buaug  = [Bu , BuEnr];
                Buaugv  = [Bu , Gkson];
        
                % Numerical integration of the sub-matrices
                Kc = Kc + Buaugv' * Du * Buaug * c;

                % Numerical integration of the internal force vector
                fi = fi + Buaugv' * stress * c;

            end

            % Number of pressure dofs associated with the continuum
            nglpc = this.nglptot - this.nglpenrPf;

            % Vector of the nodal pore-pressure dofs
            if this.staticCondensation.DisplJump
                pc = this.ue(1+this.nglu:this.nglu+nglpc);
            else
                pc = this.ue(1+this.nglutot:this.nglutot+nglpc);
            end      
            pf = this.getDiscontinuityMidPlanePressure(pc);

            % Get the discontinuity stiffness matrix and internal force
            % vector
            [Kf, ~, ~, ~, ~, ~, ~, ~, fif] = this.fractureData(UeEnr(1:this.ngluenr),pf);

            % Add contribution of the fracture
            uenrdof = this.nglu+1:this.nglutot;
            Kc(uenrdof,uenrdof) = Kc(uenrdof,uenrdof) + Kf;
            fi(uenrdof) = fi(uenrdof) + fif;
            K = Kc;

            % % Smoothed condensed stiffness matrix
            % K = Kdd - (this.getDomainArea()/At) * Kdwi * (Kwiwi\Kdwi'); 
            % for n = 1:nNeighbors
            %     Kwnd = this.neighbors(n).Kwd;
            %     An   = this.neighbors(n).getDomainArea();
            %     cols = (n-1)*nglw_neighbors+1:n*nglw_neighbors;
            %     Kn   = -(An/At) * Kdwn(:,cols) * Kwnd;
            %     K    = [K, Kn];
            % end

            % Rotation matrix from xy to mn
            % cs = this.fracture{1}(1).m(1);
            % sn = this.fracture{1}(1).m(2); 
            % R  = [ -2.0*cs*sn , 2.0*cs*sn, cs*cs-sn*sn;
            %         sn*sn     , cs*cs    , -cs*sn    ];
            % 
            % % Estimate the smoothed displacement jump 
            % smoothedW = -R * enrStrain;
            % 
            % % Update the displacement jump and traction at the fracture
            % this.fracture{1}(1).updateTractionStress(smoothedW);

        end
        
        % Assemble the element matrices and vectors
        function [Ke, Ce, fe] = assembleElementData(this, K, H, S, ~, ~, fi, fp)

            % Initialize element internal force vector
            fe = zeros(this.ngle, 1);

            % Initialize auxiliary zero matrices
            Ouu = zeros(this.nglu, this.nglu);
            Opu = zeros(this.nglptot, this.nglu);

            % Local degrees of freedom
            uDofs = 1:this.nglu;
            wDofs = 1+this.nglu:this.nglutot;

            % Condensed stiffness matrix
            Kc = K(uDofs,uDofs) - K(uDofs,wDofs) * (K(wDofs,wDofs) \ K(wDofs,uDofs));

            % Assemble internal force vector
            fe(uDofs) = fi(uDofs);
            fe(1+length(uDofs):end) = fp;
            

            % Assemble element matrices
            Ke = [  Kc , Opu';
                   Opu ,  H ];

            Ce = [ Ouu ,  Opu';
                   Opu ,   S ];

        end

        %------------------------------------------------------------------
        function Benr = getBEnrMatrix(this)

            % Compute the B matrix at the int. point and the detJ
            [dNdx, ~] = this.shape.dNdxMatrix(this.node,this.intPoint(1).X);

            % Compute the enriched B-matrix
            [~, Benr] = this.enhancedCompatibilityMtrx(dNdx,this.intPoint(1).X);

        end

        % -----------------------------------------------------------------
        % Compute the enriched B-matrix
        function [BpEnr,BuEnr] = enhancedCompatibilityMtrx(this, B, ~)

            % Initialize the enhanced shape function matrix
            BpEnr = zeros(size(B,1),this.nglpenrDp);
            BuEnr = zeros(size(this.shape.BMatrix(B),1),this.ngluenr);

            for i = 1:this.nfrac
    
                % Compute the Heaviside matrix
                Hd = this.heavisideMtrx(this.fracture{i});

                % Dofs for assemblage
                np = this.nglpenrDp/this.nfrac;
                nu = this.ngluenr/this.nfrac;
                dofp_i = np*(i - 1) + 1;
                dofp_f = np * i;
                dofu_i = nu*(i - 1) + 1;
                dofu_f = nu * i;

                % Enriched B
                Benr = -B*Hd;
    
                % Compute the enriched B-matrix
                BpEnr(:,dofp_i:dofp_f) = Benr * this.Mp{i};
                BuEnr(:,dofu_i:dofu_f) = this.shape.BMatrix(Benr) * this.Mu{i};

            end

        end

        % % -----------------------------------------------------------------
        % % Influence area factor
        % function f = areaFactor(this)
        % 
        % 
        % end
        

    end

end