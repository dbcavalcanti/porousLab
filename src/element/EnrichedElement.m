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
classdef EnrichedElement < RegularElement
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        idEnr              = [];            % Vector identifying the intersections
        fracture           = [];            % Object fracture
        nfrac              = 0;             % Number of main fractures
        gluenr             = [];            % Vector of the enhancement degrees of freedom
        glpenr             = [];            % Vector of the enhancement degrees of freedom
        glpenrDp           = [];            % Vector of the enhancement degrees of freedom
        glpenrPf           = [];            % Vector of the enhancement degrees of freedom
        ngluenr            = 0;             % Number of enhancement dof
        ngluenrFrac        = [];            % Number of jump pressure dofs per main fracture
        nglpenr            = 0;             % Number of enhancement dof
        nglutot            = 0;             % Number of enhancement dof
        nglptot            = 0;             % Number of enhancement dof
        nglpenrDp          = 0;             % Number of jump pressure dofs
        nglpenrDpFrac      = [];            % Number of jump pressure dofs per main fracture
        nglpenrPf          = 0;             % Number of longitudinal pressure dofs
        enrVar             = 'w';           % Enrichment displacement variable: 'w' or 'alpha'
        stretch            = false;         % Flag to indicate if the stretch part of the displacement mapping matrix will be considered
        subDivInt          = false;         % Flag to apply a subdivision of the domain to define the quadrature points
        jumpOrder          = 1;             % Order of the interpolation of the jump displacement field
        staticCondensation = false          % Flag to apply a static condensation of the additional dofs
        CoeffStress        = [];            % Coefficient of the polynomial approximation of the projection matrix
        fractIntersectionFlag = [];
        fractIntersectionId   = [];
        nFractIntersectionDof = 0;
        totNumberFracSegments = 0;
        EAS = false;
        Mp = [];
        Mu = [];
    end
    %% Constructor method
    methods
        function this = EnrichedElement(type, node, elem, anm, t, matModel, mat, intOrder, glu, glp, fracture, subDivInt, stretch, enrVar, jumpOrder, staticCondensation, fractIntersectionFlag, nFractIntersectionDof, fractIntersectionId, EAS)
            this = this@RegularElement(type, node, elem, anm, t, matModel, mat, intOrder, glu, glp);
            if nargin > 6
                this.isEnriched            = true;         
                this.fracture              = fracture;
                this.nfrac                 = size(fracture,1);
                this.fractIntersectionFlag = fractIntersectionFlag;
                this.fractIntersectionId   = fractIntersectionId;
                this.nFractIntersectionDof = nFractIntersectionDof;
                this.subDivInt             = subDivInt;
                this.stretch               = stretch;
                this.enrVar                = enrVar;
                this.jumpOrder             = jumpOrder;
                this.staticCondensation    = staticCondensation;
                this.EAS                   = EAS; 
                this.initializeEnrichedDofData();
            end
        end
    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        function initializeEnrichedDofData(this)

            % Initialize vectors with enriched dofs
            this.gluenr   = [];
            this.glpenrDp = [];
            this.glpenrPf = [];

            % Initialize vector with the number of jump dofs in each (main) fracture
            this.ngluenrFrac   = zeros(this.nfrac,1);
            this.nglpenrDpFrac = zeros(this.nfrac,1);

            % Get the enrichment dofs doing a loop through the fractures
            for i = 1:this.nfrac
                nfracSeg = length(this.fracture{i});
                this.totNumberFracSegments = this.totNumberFracSegments + nfracSeg;
                countDp = 0;
                countDu = 0;
                for k = 1:nfracSeg
                    this.gluenr   = [this.gluenr, this.fracture{i}(k).glu];
                    this.glpenrDp = [this.glpenrDp, this.fracture{i}(k).glDp];
                    this.glpenrPf = [this.glpenrPf, this.fracture{i}(k).glpf];
                    countDp = countDp + length(this.fracture{i}(k).glDp);
                    countDu = countDu + length(this.fracture{i}(k).glu);
                end   
                this.nglpenrDpFrac(i) = countDp;
                this.ngluenrFrac(i)   = countDu;
            end

            % Concatenate the both types of enrichment associated with the
            % pressure field
            this.glpenr    = [this.glpenrDp , this.glpenrPf]; 

            % Compute the number of enrichment dofs per type
            if this.staticCondensation.DisplJump
                this.ngluenr = this.nfrac * this.fracture{1}(1).nglu;
            else
                this.ngluenr   = length(this.gluenr);
            end
            this.nglpenrDp = length(this.glpenrDp);
            this.nglpenrPf = length(this.glpenrPf);
            this.nglpenr   = length(this.glpenr);

            % Update total number of dofs
            this.nglutot   = this.nglu + this.ngluenr;
            this.nglptot   = this.nglp + this.nglpenr;

        end

        %------------------------------------------------------------------
        function fillDofData(this)

            % Assemble the element dofs
            this.gle       = [this.glu , this.gluenr, this.glp , this.glpenrDp , this.glpenrPf];   
            this.ngle      = length(this.gle);
            
        end

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(this.intOrder,this);

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                if strcmp(this.matModel,'poroElastic')
                    constModel = MaterialHydroMechanical_Elastic(this.mat, this.anm);
                end
                intPts(i) = IntPoint(X(:,i),w(i),this.anm, constModel);
            end
            this.intPoint = intPts;

            this.CoeffStress = this.getPolynomialCoeffs();

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
            
            % Compute element matrices and internal force vector
            if this.staticCondensation.DisplJump
                [~, K, H, S, Q1, Q2, fi, fp] = this.solveLocalEquilibrium();
            else
                % Get element displacement jumps
                UeEnr = this.getDisplacementJumps();

                % Fill element matrices and internal force vector
                [K, H, S, Q1, Q2, fi, fp] = this.elementEnrichmentData(UeEnr);
            end

            % Assemble the element matrices
            [Ke, Ce, fe] = this.assembleElementData(K, H, S, Q1, Q2, fi, fp);

        end
        
        %------------------------------------------------------------------
        % Assemble the element matrices and vectors
        function [Ke, Ce, fe] = assembleElementData(this, K, H, S, Q1, Q2, fi, fp)

            % Initialize element internal force vector
            fe = zeros(this.ngle, 1);

            if this.staticCondensation.DisplJump

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

            else

                % Initialize auxiliary zero matrices
                Ouu = zeros(this.nglutot, this.nglutot);
                Opu = zeros(this.nglptot, this.nglutot);

                % Stiffness matrix
                Kc = K;

                % Assemble internal force vector
                fe(1:this.nglutot)     = fi;
                fe(1+this.nglutot:end) = fp;

            end
            

            % Assemble element matrices
%             Ke = [  Kc , -Q1;
%                    Opu ,  H ];
%             
%             Ce = [ Ouu ,  Opu';
%                     Q2' ,   S ];

            Ke = [  Kc , Opu';
                   Opu ,  H ];

            Ce = [ Ouu ,  Opu';
                   Opu ,   S ];

        end


        %------------------------------------------------------------------
        function [K, H, S, Q1, Q2, fi, fp] = elementEnrichmentData(this, UeEnr)

            % Number of pressure dofs associated with the continuum
            nglpc = this.nglptot - this.nglpenrPf;

            % Initialize the sub-matrices of the continuum
            Kc = zeros(this.nglutot, this.nglutot);  % Stiffness matrix 
            Hc = zeros(nglpc       , nglpc);         % Fluid-flow matrix 
            Sc = zeros(nglpc       , nglpc);         % Compressibility matrix
            Qc = zeros(this.nglutot, nglpc);         % HM coupling matrix

            % Initialize the internal force sub-vectors
            fic = zeros(this.nglutot, 1);

            % Initialize 2D identity vector
            m = [1.0 ; 1.0 ; 0.0];

            % Vector of the nodal pore-pressure dofs
            if this.staticCondensation.DisplJump
                pc = this.ue(1+this.nglu:this.nglu+nglpc);
            else
                pc = this.ue(1+this.nglutot:this.nglutot+nglpc);
            end      
            pf = this.getDiscontinuityMidPlanePressure(pc);
            P  = [pc; pf];

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);
                
                % Enhanced shape function matrix
                NpEnr = this.enhancedShapeFncMtrx(Np,this.intPoint(i).X);

                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.shape.BMatrix(Bp);

                % Compute the enriched B-matrix
                [BpEnr, BuEnr] = this.enhancedCompatibilityMtrx(Bp,this.intPoint(i).X);

                % Compute enhanced strain discretization matrix
                G = this.enhancedStaticMatrix(this.fracture{1}(1),this.intPoint(i).X);
                
                % Compute the increment of the strain vector
                strain = Bu*this.ue(1:this.nglu) + BuEnr*UeEnr;

                % Compute the stress vector and the constitutive matrix
                [stress,Du] = this.intPoint(i).constitutiveModel(strain);

                % Compute the stress vector and the constitutive matrix
                Dh = this.intPoint(i).constitutiveMdl.permeabilityMtrx();

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibilityCoeff();

                % Get Biot's coefficient
                biot = this.intPoint(i).constitutiveMdl.biotCoeff();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;

                % Construct augmented matrices
                Buaug  = [Bu , BuEnr];
                Buaugv = [Bu , G];
                Npaug  = [Np , NpEnr];
                Bpaug  = [Bp , BpEnr];
        
                % Numerical integration of the sub-matrices
%                 Kc = Kc + Buaug' * Du   * Buaug * c;
                Kc = Kc + Buaugv' * Du   * Buaug * c;
                Hc = Hc + Bpaug' * Dh   * Bpaug * c;
                Sc = Sc + Npaug' * comp * Npaug * c;
                Qc = Qc + Buaugv' * biot *   m   * Npaug * c;

                % Pore pressure at the integration point
                pIP = Npaug * pc;

                % Numerical integration of the internal force vector
%                 fic = fic + Buaug' * (stress - biot * pIP * m) * c;
                fic = fic + Buaugv' * (stress - biot * pIP * m) * c;

            end

            % Get the discontinuity stiffness matrix and internal force
            % vector
            [Kf, L1, L2, L3, Hf, Qf, Sf, Ef, fif, Ks] = this.fractureData(UeEnr,pf);

            % Add contribution of the fracture to the system
            [K, H, S, Q1, Q2, fi] = this.addFractureContribution(Kc, Hc, Sc, Qc, Kf, Hf, Qf, Sf, Ef, L1, L2, L3, fic, fif, Ks, pc);

            % Hydraulic part of the internal force vector
            if this.staticCondensation.Pf == true
                fp = H*pc;
            else
                fp = H*P;
            end
        end

        %------------------------------------------------------------------
        function [UeEnr, K, H, S, Q1, Q2, fi, fp] = solveLocalEquilibrium(this)

            % Parameters for the iterative process
            tol   = 1.0e-5;
            maxIt = 20;
            iter  = 1;

            % Number of displacement jump dofs
            UeEnr = this.getNodalDisplacementJumps();

            % Discontinuity local dofs
            dDofs = this.nglu+1:this.nglutot;

            % Solve the local equilibrium equation
            while (iter < maxIt)

                % Evaluate the system
                [K, H, S, Q1, Q2, fi, fp] = this.elementEnrichmentData(UeEnr);

                % Check convergence
                if norm(fi(dDofs)) < tol
                    break
                end

                % Update the displacement jumps
                dUeEnr = -K(dDofs,dDofs)\fi(dDofs);
                UeEnr  = UeEnr + dUeEnr;

                % Update the iterative counter
                iter = iter + 1;

            end

        end

        % -----------------------------------------------------------------
        function UeEnr = getNodalDisplacementJumps(this)

            % Number of displacement jump dofs
            UeEnr = zeros(this.ngluenr,1);

            count = 1;
            for i = 1:this.nfrac
                % Number of fracture segments
                nfracSeg = length(this.fracture{i});
                for j = 1:nfracSeg
                    R = this.fracture{i}(j).rotationPointMtrx();
                    for k = 1:this.fracture{i}(j).nIntPoints
                        ul = this.fracture{i}(j).intPoint(i).strainOld;
                        UeEnr(count:count+this.fracture{i}(j).ndof_nd-1) = R'*ul;
                        count = count + this.fracture{i}(j).ndof_nd;
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        % Function to get the matrices and vectors of the fractures
        function [Kf, Lm, Lmf, Lf, Hf, Qf, Sf, Ef, fif, Ks] = fractureData(this,UeEnr,pf)

            % Number of pressure dofs associated with the continuum
            nglpc = this.nglptot - this.nglpenrPf;

            % Number of dofs per fracture
            countU  = 1;
            countP  = 1;
            countPf = 1;

            % Initialize the sub-matrices of the continuum
            Kf  = zeros(this.ngluenr,this.ngluenr);           % Stiffness matrix
            Ks  = zeros(this.ngluenr,this.nglutot);           % Stabilization matrix
            Hf  = zeros(length(pf),  length(pf));             % Fluid-flow matrix 
            Sf  = zeros(length(pf),  length(pf));             % Compressibility matrix
            Lm  = zeros(nglpc,       nglpc);                  % Fluid-flow coupling matrix
            Lmf = zeros(nglpc,       length(pf));             % Fluid-flow coupling matrix
            Lf  = zeros(length(pf),  length(pf));             % Fluid-flow coupling matrix
            Qf  = zeros(this.ngluenr,length(pf));             % HM coupling matrix
            Ef  = zeros(this.nglutot,length(pf));             % HM coupling matrix
            fif = zeros(this.ngluenr,1);                      % Internal force vector

            for i = 1:this.nfrac
                
                % Number of fracture segments
                nfracSeg = length(this.fracture{i});

                for j = 1:nfracSeg

                    % Get the dofs of segment j of fracture i (in the element local numbering)
                    udofF  = countU:countU+this.fracture{i}(j).nglu-1;
                    pdofF  = countP:countP+this.fracture{i}(j).ndofDp-1;
                    pfdofF = countPf:countPf+this.fracture{i}(j).ndofPf-1;

                    % Get the matrices of segment j of fracture i
                    [Kf_i, L1_i, L2_i, L3_i, Hf_i, Qf_i, Sf_i, Ef_i, fif_i, Ks_i] = this.fracture{i}(j).elementData(UeEnr(udofF),pf(pfdofF),i,this);

                    % Assemble
                    Kf(udofF,udofF)    = Kf(udofF,udofF)    + Kf_i;   
                    Ks(udofF,:)        = Ks(udofF,:)        + Ks_i; 
                    Hf(pfdofF,pfdofF)  = Hf(pfdofF,pfdofF)  + Hf_i;
                    Sf(pfdofF,pfdofF)  = Sf(pfdofF,pfdofF)  + Sf_i;
                    Lm = Lm  + L1_i;
                    Lmf(:,pfdofF) = Lmf(:,pfdofF) + L2_i;
                    Lf(pfdofF,pfdofF)  = Lf(pfdofF,pfdofF)  + L3_i;
                    Qf(udofF,:)   = Qf(udofF,:)   + Qf_i;
                    fif(udofF,1)       = fif(udofF,1)       + fif_i;
                    % Ef  = 

                    % Update counters
                    countU  = countU + this.fracture{i}(j).nglu;
                    countP  = countP + this.fracture{i}(j).ndofDp;
                    countPf = countPf + this.fracture{i}(j).ndofPf;

                end  
            end
        end

        % -----------------------------------------------------------------
        % Function to assemble the element stiffness matrix and internal
        % force vector.
        % The assembly is based on the flag to apply a static condensation
        % or not.
        function [Ke, He, Se, Qe1, Qe2, fi] = addFractureContribution(this, Ke, Hc, Sc, Qc, Kf, Hf, Qf, Sf, Ef, L1, L2, L3, fi, fif, Ks, pc)

            % -- Mechanical internal force vector -------------------------
            uenrdof = this.nglu+1:this.nglutot;
            fi(uenrdof) = fi(uenrdof) + fif;

            % -- Stiffness matrix -----------------------------------------
            Ke(uenrdof,uenrdof) = Ke(uenrdof,uenrdof) + Kf;

            % -- Matrices of the continuity equations ---------------------
            % Add the transversal flow coupling contribution
            Hc = Hc + L1;    % to the porous-media matrix
            Hf = Hf + L3;    % to the discontinuity matrix

            Quf = [zeros(this.nglu,2); Qf]; %Valid for a single fracture

            if this.staticCondensation.Pf == false
               
                % Fluid-flow matrix
                He = [  Hc , -L2;
                       -L2',  Hf];
    
                % -- Compressibiliy matrix --------------------------------
                Se = zeros(size(He));
                
                % Add contribution of the porous media
                dof = 1:(this.nglp+this.nglpenrDp);
                Se(dof,dof) = Sc;
                
                % Add contribution of the fracture
                dofTot = 1:size(He,1);
                dof = setdiff(dofTot,dof);
                Se(dof,dof) = Sf;
    
                % -- HM-Coupling matrix -----------------------------------
                % Qe1 = zeros(this.nglutot,this.nglptot);
                % Qe2 = zeros(this.nglutot,this.nglptot);
                Qe1  = [Qc , Quf];
                % Qe2  = [Qc' ; (Quf' + Ef)];
                Qe2  = Qe1;

            % Static condensation of the longitudinal pressure dofs -------
            elseif this.staticCondensation.Pf == true
            
                % Compute the condensation matrix
                R = this.computeSystemCondensationMtrxPf(pc);

                % Apply the static condensation
                He  = Hc + R'*Hf*R - L2*R - R'*L2';
                Se  = Sc + R'*Sf*R;
                Qe1 = Qc + Quf*R;
                % Qe2 = Qc' + R'*(Ef + Quf');
                Qe2 = Qe1;
                % Qe1 = zeros(this.nglutot,this.nglptot);
                % Qe2 = zeros(this.nglutot,this.nglptot);

            end

        end

        %------------------------------------------------------------------
        % Function to initialize the stress
        function initialGeostaticStresses(this, K0, yTop, patm)

            % Initialize the stresses at the continuum element IP
            initialGeostaticStresses@RegularElement(this, K0, yTop, patm);

            % Initialize the stresses at the fracture IP
            for i = 1:this.nfrac
                % Get the number of fracture sub-segments
                nfracSeg = length(this.fracture{i});
                for j = 1:nfracSeg
                    this.fracture{i}(j).initialGeostaticStresses(this, K0, yTop, patm);
                end
            end

        end

        % -----------------------------------------------------------------
        % Function to compute the static condensation matrix of Pf
        function R = computeSystemCondensationMtrxPf(this,pc)

            % Number of pressure dofs associated with the continuum
            nglpc = this.nglptot - this.nglpenrPf;

            % Initialize the condensation matrix
            R = zeros(2*this.totNumberFracSegments,nglpc);

            % Fill condensation matrix
            count = 1;
            for i = 1:this.nfrac

                % Get the number of fracture sub-segments
                nfracSeg = length(this.fracture{i});

                for j = 1:nfracSeg      

                    for k = 1:this.fracture{i}(j).nnd_el

                        % Get the node of the discontinuity
                        Xf = this.fracture{i}(j).node(k,:);
    
                        % Get the natural coordinates of the node
                        Xn = this.shape.coordCartesianToNatural(this.node,Xf);
    
                        % Regular shape function matrix of the continuum
                        N = this.shape.shapeFncMtrx(Xn);
    
                        % Get the ids of the fracture connected to this
                        % node
                        connectedFract = this.getConnectedFracturesToNode(i,j,k);

                        if isfield(this.staticCondensation,{'PfCriteria'})==false
                            Nenr = this.middleEnhancedShapeFncMtrx(N,Xn,connectedFract);
                        elseif strcmp(this.staticCondensation.PfCriteria,'mean')
                            Nenr = this.middleEnhancedShapeFncMtrx(N,Xn,connectedFract);
                        elseif strcmp(this.staticCondensation.PfCriteria,'max')
                            Ntop = this.topEnhancedShapeFncMtrx(N,Xn,i);
                            Nbot = this.bottomEnhancedShapeFncMtrx(N,Xn,i);
                            pt = [N, Ntop]*pc;
                            pb = [N, Nbot]*pc;
                            if pt > pb
                                Nenr = Ntop;
                            elseif pt < pb
                                Nenr = Nbot;
                            elseif pt == pb
                                Nenr = 0.5*(Nbot + Ntop);
                            end
                        end
                    
                        % Assemble the condensation of the dof
                        R(count,:) = [N , Nenr];

                        count = count + 1;
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        function connectedFract = getConnectedFracturesToNode(this,mainFractId,fractSegId,nodeFractSeg)

            % The main fracture will always be connected to its given node
            connectedFract = mainFractId;

            % Check if this node corresponds to an intersection with
            % another fracture
            if (this.fractIntersectionFlag==true)

                % Get the nodes of the main fracture that corresponds to a
                % intersection
                intersectionNodes = this.fractIntersectionId(mainFractId,:);
                [~,intIdFracj,ndIntFracti] = find(intersectionNodes);

                % In case it corresponds, add it to the array
                if (this.fracture{mainFractId}(fractSegId).connect(nodeFractSeg) == ndIntFracti) %|| (this.fracture{mainFractId}(fractSegId).connect(nodeFractSeg) == ndIntFracti+1)
                    connectedFract = [mainFractId,intIdFracj];
                end
            end
        end

        % -----------------------------------------------------------------
        % Function to update the state variables
        function updateStateVar(this)

            % Update the state variables at the continuum element IPnts
            updateStateVar@RegularElement(this);

            % Update the state variables at the fracture element IPnts
            for i = 1:this.nfrac
                for j = 1:length(this.fracture{i})
                    this.fracture{i}(j).updateStateVar();
                end
            end

        end

        % -----------------------------------------------------------------
        % Function to compute the enrichment dofs.
        % When a static condensation of the additional dofs is applied,
        % this increments must be computed.
        function UeEnr = getDisplacementJumps(this)

            UeEnr = this.ue((1+this.nglu):((this.nglu+this.ngluenr)));

        end

        % -----------------------------------------------------------------
        % Function to compute the increment of the enrichment dofs.
        % When a static condensation of the additional dofs is applied,
        % this increments must be computed.
        function pf = getDiscontinuityMidPlanePressure(this,pc)
            
            if this.staticCondensation.Pf == true

                % Initialize pressure vector
                pf = zeros(2*this.totNumberFracSegments,1);

                count = 1;
                for i = 1:this.nfrac

                    % Get the number of fracture sub-segments
                    nfracSeg = length(this.fracture{i});

                    for j = 1:nfracSeg

                        % Get fracture nodes
                        for k = 1:this.fracture{i}(j).nnd_el
        
                            % Get the node of the discontinuity
                            Xf = this.fracture{i}(j).node(k,:);
        
                            % Get the natural coordinates of the node
                            Xn = this.shape.coordCartesianToNatural(this.node,Xf);
        
                            % Regular shape function matrix of the continuum
                            N = this.shape.shapeFncMtrx(Xn);
                            
                            % Get the ids of the fracture connected to this
                            % node
                            connectedFract = this.getConnectedFracturesToNode(i,j,k);
                            
                            if isfield(this.staticCondensation,{'PfCriteria'})==false
                                Nenr = this.middleEnhancedShapeFncMtrx(N,Xn,connectedFract);
                            elseif strcmp(this.staticCondensation.PfCriteria,'mean')
                                Nenr = this.middleEnhancedShapeFncMtrx(N,Xn,connectedFract);
                            elseif strcmp(this.staticCondensation.PfCriteria,'max')
                                Ntop = this.topEnhancedShapeFncMtrx(N,Xn,i);
                                Nbot = this.bottomEnhancedShapeFncMtrx(N,Xn,i);
                                pt = [N, Ntop]*pc;
                                pb = [N, Nbot]*pc;
                                if pt > pb
                                    Nenr = Ntop;
                                elseif pt < pb
                                    Nenr = Nbot;
                                elseif pt == pb
                                    Nenr = 0.5*(Nbot + Ntop);
                                end
                            end

                            % Compute the mean pressure at the discontinuity node
                            pf(count) = [N , Nenr] * pc;
    
                            count = count + 1;
                        end
                    end
                end

            else
                
                % Number of pressure dofs associated with the continuum
                nglpc = length(pc);

                if this.staticCondensation.DisplJump
                    pf = this.ue(1+this.nglu+nglpc:end);
                else
                    pf = this.ue(1+this.nglutot+nglpc:end);
                end  
                
            end

        end

        % -----------------------------------------------------------------
        % Compute the enhanced shape function matrix 
        function Nenr = enhancedShapeFncMtrx(this, N, Xn)

            % Integration point in the cartesian coordinates
            X = this.shape.coordNaturalToCartesian(this.node, Xn);
            
            % Initialize the enhanced shape function matrix
            Nenr = zeros(size(N,1),this.nglpenrDp);
            
            for i = 1:this.nfrac

                % Evaluate the Heaviside function at the point X
                h = this.heavisideFnc(this.fracture{i},X);
    
                % Compute the Heaviside matrix
                Hd = this.heavisideMtrx(this.fracture{i});

                % Identity matrix
                I = eye(size(Hd,1));

                % Dofs for assemblage
                dof_i = this.nglpenrDpFrac(i)*(i - 1) + 1;
                dof_f = this.nglpenrDpFrac(i) * i;
    
                % Compute the matrix Gr
                Nenr(:,dof_i:dof_f) = N*(h*I - Hd)*this.Mp{i};

            end
        end

        % -----------------------------------------------------------------
        % Compute the enhanced shape function matrix 
        function Nenr = displacementEnhancedShapeFncMtrx(this, N, Xn)

            % Integration point in the cartesian coordinates
            X = this.shape.coordNaturalToCartesian(this.node, Xn);
            
            % Initialize the enhanced shape function matrix
            Nenr = zeros(2,this.ngluenr);
            
            for i = 1:this.nfrac

                % Evaluate the Heaviside function at the point X
                h = this.heavisideFnc(this.fracture{i},X);
    
                % Compute the Heaviside matrix
                Hd = this.heavisideMtrx(this.fracture{i});

                % Identity matrix
                I = eye(size(Hd,1));

                % Assemble the element mapping matrix 
                if this.staticCondensation.DisplJump
                    dof_i = this.fracture{1}(1).nglu*(i - 1) + 1;
                    dof_f = this.fracture{1}(1).nglu * i;
                else
                    dof_i = this.ngluenrFrac(i)*(i - 1) + 1;
                    dof_f = this.ngluenrFrac(i) * i;
                end
                
    
                % Compute the matrix Gr
                Nenr(:,dof_i:dof_f) = this.shape.NuMtrx(N*(h*I - Hd))*this.Mu{i};

            end
        end

        % -----------------------------------------------------------------
        % Compute the enhanced shape function matrix at the region Omega^+
        function Nenr = topEnhancedShapeFncMtrx(this, N, Xn, idFrac)

            % Integration point in the cartesian coordinates
            X = this.shape.coordNaturalToCartesian(this.node, Xn);
            
            % Initialize the enhanced shape function matrix
            Nenr = zeros(size(N,1),this.nglpenrDp);
            
            for i = 1:this.nfrac
                if i == idFrac
                    h = 1.0;
                else
                    % Evaluate the Heaviside function at the point X
                    h = this.heavisideFnc(this.fracture{i},X);
                end
    
                % Compute the Heaviside matrix
                Hd = this.heavisideMtrx(this.fracture{i});

                % Identity matrix
                I = eye(size(Hd,1));

                % Dofs for assemblage
                dof_i = this.nglpenrDpFrac(i)*(i - 1) + 1;
                dof_f = this.nglpenrDpFrac(i) * i;
    
                % Compute the matrix Gr
                Nenr(:,dof_i:dof_f) = N*(h*I - Hd)*this.Mp{i};

            end
        end

        % -----------------------------------------------------------------
        % Compute the enhanced shape function matrix at the region Omega^-
        function Nenr = bottomEnhancedShapeFncMtrx(this, N, Xn, idFrac)

            % Integration point in the cartesian coordinates
            X = this.shape.coordNaturalToCartesian(this.node, Xn);
            
            % Initialize the enhanced shape function matrix
            Nenr = zeros(size(N,1),this.nglpenrDp);
            
            for i = 1:this.nfrac
                if i == idFrac
                    h = 0.0;
                else
                    % Evaluate the Heaviside function at the point X
                    h = this.heavisideFnc(this.fracture{i},X);
                end
    
                % Compute the Heaviside matrix
                Hd = this.heavisideMtrx(this.fracture{i});

                % Identity matrix
                I = eye(size(Hd,1));

                % Dofs for assemblage
                dof_i = this.nglpenrDpFrac(i)*(i - 1) + 1;
                dof_f = this.nglpenrDpFrac(i) * i;
    
                % Compute the matrix Gr
                Nenr(:,dof_i:dof_f) = N*(h*I - Hd)*this.Mp{i};

            end
        end

        % -----------------------------------------------------------------
        % Compute the mean enhanced shape function matrix 
        function Nenr = middleEnhancedShapeFncMtrx(this,N, Xn, connectedFrac)

            % Integration point in the cartesian coordinates
            X = this.shape.coordNaturalToCartesian(this.node, Xn);
            
            % Initialize the enhanced shape function matrix
            Nenr = zeros(size(N,1),this.nglpenrDp);
            
            % Fill the enriched shape function matrix
            for i = 1:this.nfrac

                if ismember(i,connectedFrac)
                    h = 0.5;
                else
                    % Evaluate the Heaviside function at the point X
                    h = this.heavisideFnc(this.fracture{i}(1),X);
                end

                % Compute the Heaviside matrix
                Hd = this.heavisideMtrx(this.fracture{i}(1));

                % Identity matrix
                I = eye(size(Hd,1));

                % Dofs for assemblage
                dof_i = this.nglpenrDpFrac(i)*(i - 1) + 1;
                dof_f = this.nglpenrDpFrac(i) * i;
    
                % Compute the matrix Gr
                Nenr(:,dof_i:dof_f) = N*(h*I - Hd)*this.Mp{i};

            end
        end

        % -----------------------------------------------------------------
        % Compute the mean enhanced shape function matrix 
        function Nstab = StabEnhancedShapeFncMtrx(this,idfrac)

            % Initialize the enhanced shape function matrix
            Ae = 0.0;
            Nstab = zeros(2,this.ngluenr);

            for i = 1:this.nIntPoints

                [~, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;

                % Shape function matrix of the continuum
                N = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Integrate
                Nstab = Nstab + this.shape.NuMtrx(N) * this.Mu{idfrac} * c;
                Ae    = Ae + c;

            end
            Nstab = Nstab / Ae;
        end

        % -----------------------------------------------------------------
        % Compute the enriched B-matrix
        function [BpEnr,BuEnr] = enhancedCompatibilityMtrx(this, B, Xn)

            % Integration point in the cartesian coordinates
            X = this.shape.coordNaturalToCartesian(this.node, Xn);

            % Initialize the enhanced shape function matrix
            BpEnr = zeros(size(B,1),this.nglpenrDp);
            BuEnr = zeros(size(this.shape.BMatrix(B),1),this.ngluenr);

            for i = 1:this.nfrac

                % Evaluate the Heaviside function at the point X
                h = this.heavisideFnc(this.fracture{i},X);
    
                % Compute the Heaviside matrix
                Hd = this.heavisideMtrx(this.fracture{i});

                % Identity matrix
                I = eye(size(Hd,1));

                % Mapping matrices

                % Dofs for assemblage
                np = this.nglpenrDp/this.nfrac;
                nu = this.ngluenr/this.nfrac;
                dofp_i = np*(i - 1) + 1;
                dofp_f = np * i;
                dofu_i = nu*(i - 1) + 1;
                dofu_f = nu * i;

                % Enriched B
                Benr = B*(h*I - Hd);
    
                % Compute the enriched B-matrix
                BpEnr(:,dofp_i:dofp_f) = Benr * this.Mp{i};
                BuEnr(:,dofu_i:dofu_f) = this.shape.BMatrix(Benr) * this.Mu{i};

            end

        end

        % -----------------------------------------------------------------
        % Compute the mean enriched B-matrix: (Btop + Bbot)*0.5
        function Benr = enhancedMiddleCompatibilityMtrx(this, B, ~)

            % Compute the Heaviside matrix
            Hd = this.heavisideMtrx();

            % Identity matrix
            I = eye(size(Hd,1));

            % Compute the enriched B-matrix
            Benr = B*(0.5*I - Hd);

        end

        % -----------------------------------------------------------------
        function [stress_gamma, DBaug] = weightedStress(this,idFrac,X)

            % Get the natural coordinates of the node
            Xn = this.shape.coordCartesianToNatural(this.node,X);

            % Compute the B matrix at the int. point and the detJ
            dNdX = this.shape.dNdxMatrix(this.node,Xn);

            % Assemble the B-matrix for the mechanical part
            B = this.shape.BMatrix(dNdX);

            % weighted enriched B-matrix
            Benrg = this.weightedEnhancedCompatibilityMtrx(dNdX, idFrac);

            % Augmented B-matrix
            Baug = [B , Benrg];

            % Get the elastic constitutive matrix
            D = this.intPoint(1).constitutiveMdl.elasticConstitutiveMtrx();

            % Current total displacement
            ue = this.ue(1:this.nglutot);

            % Compute the stress
            DBaug        = D * Baug;
            stress_gamma = DBaug * ue;


        end

        % -----------------------------------------------------------------
        function Benrg = weightedEnhancedCompatibilityMtrx(this, dNdX, idFrac)

            % Initialize the enhanced shape function matrix
            Benrg = zeros(size(this.shape.BMatrix(dNdX),1),this.ngluenr);

            % Compute the weight associated with the fracture idFrac
            gammaTop = this.stabilizationWeight(this.fracture{idFrac});

            for i = 1:this.nfrac

                if i == idFrac
                    htop = 1.0;
                    hbot = 0.0;
                else
                    % Evaluate the Heaviside function at the point X
                    htop = this.heavisideFnc(this.fracture{i},X);
                    hbot = htop;
                end
    
                % Compute the Heaviside matrix
                Hd = this.heavisideMtrx(this.fracture{i});

                % Identity matrix
                I = eye(size(Hd,1));

                % Dofs for assemblage
                nu = this.ngluenr/this.nfrac;
                dofu_i = nu*(i - 1) + 1;
                dofu_f = nu * i;

                % Enriched B
                Btop = dNdX*(htop * I - Hd);
                Bbot = dNdX*(hbot * I - Hd);

                % Weighted enriched matrix
                Bg = gammaTop*Btop + (1.0 - gammaTop)*Bbot;
    
                % Compute the enriched B-matrix
                Benrg(:,dofu_i:dofu_f) = this.shape.BMatrix(Bg) * this.Mu{i};

            end
        end

        % -----------------------------------------------------------------
        function G = enhancedStaticMatrix(this,fracture,Xn)
            
            % Get the fracture normal vector
            P = fracture.stressProjectionMatrix();

            % Cartesian coordinates of the given point
            X = this.shape.coordNaturalToCartesian(this.node, Xn);

            % Compute the polynomial stress approximation
            p = this.shape.polynomialStress(X);

            % Interpolation 
            gk = this.CoeffStress'*p;

            % G matrix
            G = zeros(size(P,1),length(gk)*size(P,2));
            for i = 1:length(gk)
                G(:,(i-1)*size(P,2)+1:i*size(P,2)) = -gk(i) * P;
            end
            
        end

        % -----------------------------------------------------------------
        function gamma = stabilizationWeight(this,fracture)

            % Get the area of the element domain
            A = this.getDomainArea();

            % Calculate area of the region Omega^+
            Atop = this.areaPositiveSubRegion(fracture);

            % Compute the weight associated to the sub-region Omega^+
            gamma = 0.5; %Atop/A;

        end

        % -----------------------------------------------------------------
        function Atop = areaPositiveSubRegion(this,fracture)

            % Evaluate the heaviside function at the nodes of the element
            h = this.heavisideFnc(fracture,this.node);

            % Number of fracture segments
            nfracSeg = length(fracture);

            % Get the fracture tip nodes
            fractiNodes = zeros(2);
            fractiNodes(1,:) = fracture(1).node(1,:);
            fractiNodes(2,:) = fracture(nfracSeg).node(2,:);

            % Get the nodes of the sub-region Omega^+
            Nodes = [this.node(h > 0.0,:); fractiNodes];   
            order = this.sortCounterClockWise(Nodes);

            % Calculate area of the region Omega^+
            Atop = this.calculateArea(Nodes(order,:));

        end

        % -----------------------------------------------------------------
        function computeMappingMatrices(this)
            this.Mp = cell(this.nfrac,1);
            this.Mu = cell(this.nfrac,1);
            for i = 1:this.nfrac
                % Compute the Heaviside matrix
                % Hd = this.heavisideMtrx(this.fracture{i});
                % Number of nodes on the sub-domain Omega^+
                % nNdsOplus = sum(Hd,"all");
                % Intersection mode (of the fracture with the mesh)
                % if nNdsOplus == this.nnd_el/2
                %     intMode = 'sym';
                % else
                %     intMode = 'skew';
                % end
                intMode = 'sym';
                % Store the mapping matrices
                [this.Mp{i},this.Mu{i}] = elementMappingMtrcs(this,i,intMode);
            end
        end

        % -----------------------------------------------------------------
        % Mapping matrix associated to a element. This matrix is
        % constructed by stacking by rows the mapping matrices evaluated at
        % the element's nodes.
        function [Mp,Mu] = elementMappingMtrcs(this,i,intMode)

            % Get the Poisson ratio of the material
            nu = this.mat(2);

            % Initialize the element's mapping matrix
            Mp = zeros(this.nnd_el,this.nglpenrDpFrac(i));
            if this.staticCondensation.DisplJump
                Mu = zeros(this.nnd_el*2,this.nfrac*this.fracture{1}(1).nglu);
            else
                Mu = zeros(this.nnd_el*2,this.ngluenrFrac(i));
            end

            % Get the id of the fractures j that intersects fracture i
            if isempty(this.fractIntersectionId) == false
                [~,idFractIntersect] = find(this.fractIntersectionId(i,:));
                numIntersectingFrac  = length(idFractIntersect);
            else
                idFractIntersect = [];
                numIntersectingFrac = 1;
            end

            % Construct the matrix by stacking by-rows the mapping matrix
            % evaluated at each node
            for j = 1:this.nnd_el

                % Get the coordinate of the node i
                X  = this.node(j,:);

                % Get the fracture segment that belong to the same
                % subregion of the node with coordinate X
                id = this.getFractureSegmentAtSubRegion(i,X,idFractIntersect,numIntersectingFrac);

                % Get the mapping matrix of the fracture segment
                Mpj = this.fracture{i}(id).jumpPressureTransmissionMtrx(X,intMode);
                Muj = this.fracture{i}(id).jumpDisplacementTransmissionMtrx(X,this.enrVar,this.stretch,nu,intMode);

                % Assemble
                pIni = length(this.fracture{i}(id).glDp)*(id - 1) + 1;
                pEnd = length(this.fracture{i}(id).glDp)*id;
                uIni = this.fracture{i}(id).nglu*(id - 1) + 1;
                uEnd = this.fracture{i}(id).nglu*id;

                % Assemble the element mapping matrix
                Mp(j,pIni:pEnd) = Mpj;
                Mu(2*(j-1)+1:2*j,uIni:uEnd) = Muj;
            end

        end

        % -----------------------------------------------------------------
        function idSeg = getFractureSegmentAtSubRegion(this,i,X,idFractIntersect,numIntersectingFrac)

            % Vector with the heavise function of every main
            % discontinuity evaluated at the 
            nFracSeg = length(this.fracture{i});
            
            % Identify the fracture segment that is in the same
            % subregion as the node i
            idSeg = 1;
                
            if nFracSeg > 1
                % Vector with the heavise function of every main
                % discontinuity evaluated at this node
                hX = zeros(numIntersectingFrac,1);
                for k = 1:numIntersectingFrac
                    hX(k) = this.heavisideFnc(this.fracture{idFractIntersect(k)}(1),X);
                end
                for k = 1:nFracSeg
                    hXref = zeros(numIntersectingFrac,1);
                    for ii = 1:numIntersectingFrac
                        hXref(ii) = this.heavisideFnc(this.fracture{idFractIntersect(ii)}(1),this.fracture{i}(k).Xref);
                    end
                    if norm(hXref-hX)<1.0e-15
                        idSeg = k;
                        break;
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        % Heaviside function evaluated at the point (or set of points) X,
        % wrt to the reference point of the fracture (Xref).
        function h = heavisideFnc(this,fracture,X)

            % Reference point at the fracture that crosses the element
            if this.EAS
                Xref = mean(this.node);
            else
                Xref = fracture.Xref;
            end

            % Fracture normal vector
            n = fracture.n;

            % Relative distance of the nodes of the element wrt to the reference point
            DX = X - repmat(Xref,size(X,1),1);
            
            % Heaviside function evaluated in the nodes of the element
            h = max(sign(DX*n'),0.0);
            % h = sign(DX*n');

        end

        % -----------------------------------------------------------------
        % Diagonal matrix with the Heaviside function evaluated at the
        % nodes associeted to each regular dof of the element, wrt to the
        % reference point of the fracture (Xref).
        function Hde = heavisideMtrx(this,fracture)
 
            % Heaviside function evaluated in the nodes of the element
            h = this.heavisideFnc(fracture,this.node);

            % Create the vector
            hv = repmat(h,[1,this.ndof_nd])';
            
            % Matrix with the Heaviside function evaluated in the node of the dofs
            Hde = diag(reshape(hv,numel(hv),1));

        end

        % -----------------------------------------------------------------
        % Compute the coefficients of the polynomial g^(k) that is used to 
        % approximate the submatrix matrix Gv^(k).
        % The submatrix is defined as:
        %   Gv^(k) = g^(k) * P. 
        % Where:
        %   g^(k) = c0 + c1 * xrel + c2 * yrel
        % 
        function C = getPolynomialCoeffs(this)
            
            % Gram matrix
            H = this.gramMtrx();

            % Stress interpolation vector
            S = this.fracture{1}(1).stressIntVct(this.shape,this.node);

            % Compute the coefficients
            C = H \ S;

        end

        % -----------------------------------------------------------------
        % Compute the Gramm matrix
        function H = gramMtrx(this)

            % Initialize the Gram matrix
            dim  = this.shape.getSizeGramMtrx();
            H    = zeros(dim);

            for i = 1:this.nIntPoints

                % Compute the determinant of the Jacobian matrix
                detJ = this.shape.detJacobian(this.node,this.intPoint(i).X);
        
                % Compute the integrand of the Gram Matrix
                dH = this.shape.integrandGramMtrx(this.node,this.intPoint(i).X);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
        
                % Numerical integration of the stiffness matrix Kaa
                H = H + dH * c;
            
            end

        end

        %------------------------------------------------------------------
        % Function to compute the displacement field in the element.
        function u = displacementField(this,X)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   u   : displacement vector evaluated in "X"
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);
            Nu = this.shape.NuMtrx(Nm);

            % Enriched shape function vector
            Nenr  = this.displacementEnhancedShapeFncMtrx(Nm, Xn);

            % Get enrichment dofs
            we = this.getNodalDisplacementJumps();

            % Displacement field
            regdof = 1:this.nglu;
            
            % Regular displacement field
            u = Nu*this.ue(regdof) + Nenr*we;
        
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = pressureField(this,X,ue)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   p   : pressure evaluated in "X"

            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Enriched shape function vector
            Nenr = this.enhancedShapeFncMtrx(Nm, Xn);

            % Get enrichment dofs
            Dp = this.getNodalPressureJumps();

            % Displacement field
            if this.staticCondensation.DisplJump
                regdof = (this.nglu+1):(this.nglu+this.nglp);
            else
                regdof = (this.nglutot+1):(this.nglutot+this.nglp);
            end 
            p = Nm*this.ue(regdof) + Nenr*Dp;
        
        end

        %------------------------------------------------------------------
        % Function to get the enrichment dofs in a element.
        % Considers only one discontinuity inside an element.
        function Dp = getNodalPressureJumps(this)

            if this.staticCondensation.DisplJump
                Dp = this.ue(this.nglu+this.nglp+1:this.nglu+this.nglp+this.nglpenrDp);
            else
                Dp = this.ue(this.nglutot+this.nglp+1:this.nglutot+this.nglp+this.nglpenrDp);
            end
            

        end

        %------------------------------------------------------------------
        % Compute the slip tendency parameter at a point X
        function ST = slipTendency(this,X)

            % Slip tendency at the integration points
            STip = [this.fracture{1}.intPoint(1).statevar(2); 
                    this.fracture{1}.intPoint(2).statevar(2)];

            % Normalized tangential coordinate
            s = this.fracture{1}.tangentialLocCoordinate(X);

            % Length of the discontinuity
            ld = this.fracture{1}.ld;

            % Shape function matrix
            N = [0.5-s/ld, 0.5+s/ld];
            
            % Slip tendency at X
            ST = N*STip;

        end

        %------------------------------------------------------------------
        % Compute the shear cohesive stress at X
        function ts = shearCohesiveStress(this,X)

            % Slip tendency at the integration points
            TS = [this.fracture{1}.intPoint(1).stress(1); 
                  this.fracture{1}.intPoint(2).stress(1)];

            % Normalized tangential coordinate
            s = this.fracture{1}.tangentialLocCoordinate(X);

            % Length of the discontinuity
            ld = this.fracture{1}.ld;

            % Shape function matrix
            N = [0.5-s/ld, 0.5+s/ld];
            
            % Slip tendency at X
            ts = N*TS;

        end

        %------------------------------------------------------------------
        % Compute the shear cohesive stress at X
        function tn = normalCohesiveStress(this,X)

            % Slip tendency at the integration points
            TN = [this.fracture{1}.intPoint(1).stress(2); 
                  this.fracture{1}.intPoint(2).stress(2)];

            % Normalized tangential coordinate
            s = this.fracture{1}.tangentialLocCoordinate(X);

            % Length of the discontinuity
            ld = this.fracture{1}.ld;

            % Shape function matrix
            N = [0.5-s/ld, 0.5+s/ld];
            
            % Slip tendency at X
            tn = N*TN;

        end

        %------------------------------------------------------------------
        % Function to get the nodes that will be used to plot the results.
        function [resNodes,fractNodes] = getResultNodes(this)

            resNodes = [this.node];
            
            for i = 1:this.nfrac

                % Number of fracture segments
                nfracSeg = length(this.fracture{i});

                % Get the fracture tip nodes
                fractiNodes = zeros(2);
                fractiNodes(1,:) = this.fracture{i}(1).node(1,:);
                fractiNodes(2,:) = this.fracture{i}(nfracSeg).node(2,:);

                % Order the nodes in counterclockwise order
                Nodes = [this.node; fractiNodes];
                order = this.sortCounterClockWise(Nodes);
    
                % Get the position of the fracture nodes in the order
                nf1 = find(order == size(this.node,1)+1);
                nf2 = find(order == size(this.node,1)+2);
    
                % The nodes before and after each node of the fracture define
                % the edge of the element.
                % Edge associated to the first node:
                if nf1 == 1
                    edge1 = [order(end) order(nf1+1)];
                elseif nf1 == length(order)
                    edge1 = [order(nf1-1) order(1)];
                else
                    edge1 = [order(nf1-1) order(nf1+1)];
                end
                % Edge associated to the second node:
                if nf2 == 1
                    edge2 = [order(end) order(nf2+1)];
                elseif nf2 == length(order)
                    edge2 = [order(nf2-1) order(1)];
                else
                    edge2 = [order(nf2-1) order(nf2+1)];
                 end
    
                % Length of each edge:
                len1 = sqrt((Nodes(edge1(2),1)-Nodes(edge1(1),1))^2+(Nodes(edge1(2),2)-Nodes(edge1(1),2))^2);
                len2 = sqrt((Nodes(edge2(2),1)-Nodes(edge2(1),1))^2+(Nodes(edge2(2),2)-Nodes(edge2(1),2))^2);
    
                % Compute the vectors with the origin at the nodes of the
                % fracture oriented along the edges.
                vf1 = Nodes(edge1(1),:) - Nodes(edge1(2),:);
                vf2 = Nodes(edge2(1),:) - Nodes(edge2(2),:);
                vf1 = vf1/norm(vf1);
                vf2 = vf2/norm(vf2);
    
                % Compute the increment vector:
                df1 = (len1/1000)*vf1;
                df2 = (len2/1000)*vf2;
    
                % Find a node along the Omega^plus region along each edge:
                n1_plus = fractiNodes(1,:) + df1;
                n2_plus = fractiNodes(2,:) + df2;
    
                % Find a node along the Omega^minus region along each edge:
                n1_minus = fractiNodes(1,:) - df1;
                n2_minus = fractiNodes(2,:) - df2; 
    
                % Get the order of the nodes for the element connectivity be
                % done in a counterclockwise way
                fractNodes = [n1_plus;n1_minus; n2_plus;n2_minus];
                resNodes = [resNodes; fractNodes];
                orderResNodes = this.sortCounterClockWise(resNodes);
                resNodes = resNodes(orderResNodes,:);
    
                % Nodes for the fracture result object
                orderFract = this.sortCounterClockWise(fractNodes);
                fractNodes = fractNodes(orderFract,:);
            end


        end

        %------------------------------------------------------------------
        % Function to create an array of result objects based on the
        % children elements.
        function createResults(this)
            [resNodes,fractNodes] = this.getResultNodes();
            nNodes      = size(resNodes,1);
            nFractNodes = size(fractNodes,1);
            this.result = Result(resNodes ,1:nNodes ,0.0*ones(nNodes,1) ,'');
            this.fracture{1}.result = Result(fractNodes ,1:nFractNodes ,0.0*ones(nNodes,1) ,'');
        end

    end

    
end