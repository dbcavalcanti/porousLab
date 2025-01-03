%% Fracture class
%
% This class defines a fracture element
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: December 2022
%%%
% Initially prepared for the course CIV 2801 - Fundamentos de Computacao
% Grafica, 2022, second term, Department of Civil Engineering, PUC-Rio.
%
%% Class definition
classdef Fracture < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        shape      = [];              % Object of the Shape class
        idelem     = 0;               % Identify the element that the fracture belongs
        node       = [];              % Nodes of the fem mesh
        connect    = [];              % Nodes connectivity
        m          = [];              % Tangent orientation vector
        n          = [];              % Normal orientation vector
        rs         = [];              % Rotation vector to the shear strain
        rotXY2nm   = [];              % Rotation matrix to rotate the stresses from the cartesian to the fracture local system
        ld         = 0.0;             % Fracture length
        Xref       = [];              % Reference point
        t          = 1.0;             % Thickness
        matModel   = 'elastic';       % Material model
        penal      = false;           % Flag for applying a penalization on compression
        cp         = 1.0;             % Penalization coefficient
        mat        = [];              % Vector with material properties
        nnd_el     = 2;               % Number of nodes per element
        ndof_nd    = 2;               % Number of dof per node
        ndof       = 1;               % Number of total dofs
        nglu       = 1;               % Number of displacement dofs
        nglp       = 1;               % Number of pressure dofs
        glu        = [];              % Vector of the displacement degrees of freedom
        glp        = [];              % Vector of the pressure degrees of freedom
        glDp       = [];              % Vector of the pressure degrees of freedom
        glpf       = [];              % Vector of the pressure degrees of freedom
        nIntPoints = 2;               % Number of integration points
        intPoint   = [];              % Vector with integration point objects
        result     = [];              % Result object
        intOrder   = 1;
        updatePermeability = false;
        discontinuityOpen = true;
        stabilizedCohesiveStresses = false;
        penalM = 0.0;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Fracture(nglu)
            this.nglu = nglu;
        end
    end

    methods
        function initializeFracture(this, node, elem, t, matModel, penal, cp, mat, glu, glDp, glpf, updatePermeability, discontinuityOpen, stabilizedCohesiveStresses)

            this.node     = node;
            this.connect  = elem;
            this.t        = t;
            this.matModel = matModel;
            this.penal    = penal;
            this.cp       = cp;
            this.mat      = mat;
            this.glu      = glu;
            this.glDp     = glDp;
            this.glpf     = glpf;
            this.glp      = [glDp , glpf];
            this.nglp     = length(this.glp);
            this.ndof     = this.nglu + this.nglp;
            this.shape    = Shape_Bar();
            this.updatePermeability = updatePermeability;
            this.discontinuityOpen = discontinuityOpen;
            this.stabilizedCohesiveStresses = stabilizedCohesiveStresses;

            % Initialize the geometry properties
            this.initializeGeometry();

            % Initialize the integration points
            this.initializeIntPoints();

        end
    end

    %% Abstract methods
    methods(Abstract)

        % Compute the shape function matrix
        N = shapeFncMtrx(this,xn)
        N = interpJumpShapeMtrx(this,xn, enrVar)
        
        % Compute the jump transmission matrix M
        M = jumpDisplacementTransmissionMtrx(this,X);

        % Compute the jump transmission matrix M
        M = jumpPressureTransmissionMtrx(this,X);

        % Compute the rotation matrix
        R = rotationMtrx(this);

    end  
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initialize the fracture element geometry. Computes the length,
        % the orientation vectors and the reference point.
        function initializeGeometry(this)

            % Fracture length
            dx = this.node(2,1) - this.node(1,1);
            dy = this.node(2,2) - this.node(1,2);
            this.ld = sqrt(dx.^2 + dy.^2);
            
            % Fracture orientation
            sn = dy./this.ld;
            cs = dx./this.ld;
            
            % Tangential vector to the discontinuity
            this.m = [ cs   sn];  
            
            % Normal vector to the discontinuity
            % Defined considering n = ez x m, where ez = [0 0 1]
            this.n = [-sn   cs];

            % Rotation to obtain the shear strain
            cs2 = cs*cs;
            sn2 = sn*sn;
            sncs = sn*cs;

            this.rs = [cs2 , sn2 , sncs];

            this.rotXY2nm = [ cs2 ,  sn2 ,  2.0*sncs;
                              sn2 ,  cs2 , -2.0*sncs;
                             -sncs,  sncs,  cs2-sn2];

            % Reference point
            this.Xref = 0.5*(this.node(1,:) + this.node(2,:));
            
        end

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(this.intOrder);

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints

                % Initialize the constitutive model
                if strcmp(this.matModel,'elastic')
                    constModel = MaterialHydroMechanicalInterface_Elastic(this.mat, this.penal, this.cp, this.updatePermeability, this.discontinuityOpen);
                elseif strcmp(this.matModel,'elasticRock')
                    constModel = MaterialHydroMechanicalInterface_ElasticRock(this.mat, this.penal, this.cp, this.updatePermeability, this.discontinuityOpen);
                elseif strcmp(this.matModel,'contactElasticRock')
                    constModel = MaterialHydroMechanicalInterface_ContactElasticRock(this.mat, this.penal, this.cp, this.updatePermeability, this.discontinuityOpen);
                elseif strcmp(this.matModel,'mohrCoulomb')
                    constModel = MaterialHydroMechanicalInterface_MohrCoulomb(this.mat, this.penal, this.cp, this.updatePermeability, this.discontinuityOpen);
                elseif strcmp(this.matModel,'isotropicDamageLinear')
                    constModel = MaterialHydroMechanicalInterface_IsotropicDamage_Linear(this.mat, this.penal, this.cp, this.updatePermeability, this.discontinuityOpen);
                elseif strcmp(this.matModel,'isotropicDamageExponential')
                    constModel = MaterialHydroMechanicalInterface_IsotropicDamage_Exponential(this.mat, this.penal, this.cp, this.updatePermeability, this.discontinuityOpen);
                end

                % Create the integration points
                intPts(i) = IntPoint(X(:,i),w(i),'Interface', constModel);

                % Initialize aperture
                intPts(i).statevar(1)    = constModel.getInitialAperture();
                intPts(i).statevarOld(1) = constModel.getInitialAperture();

            end

            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function updateStateVar(this)
            for i = 1:this.nIntPoints
                this.intPoint(i).updateStateVar();
                this.intPoint(i).updateStrainVct();
                this.intPoint(i).updatePlasticStrainVct();
                this.intPoint(i).updateStressVct();
            end
        end

        %------------------------------------------------------------------
        % Function to initialize the stress
        function initialGeostaticStresses(this, elem, K0, yTop, patm)

            % Loop through the integration points
            for i = 1:this.nIntPoints

                % Specific weight
                rhob = elem.intPoint(i).constitutiveMdl.rhoBulk();
                grav = elem.intPoint(i).constitutiveMdl.gravityAcc();
                gw   = elem.intPoint(i).constitutiveMdl.waterSpecificWeight();

                % Get the coordinates of the integration point
                Xn = this.intPoint(i).X;
                X  = this.shape.coordNaturalToCartesian(this.node,Xn);

                % Compute the effective stress components in the global
                % cartesian coordinate system
                sy = -patm-(yTop - X(2))*(grav*rhob - gw);
                sx = K0*sy;

                stressXY = [sx;sy;0.0];

                stressNM = this.rotXY2nm*stressXY;

                % Store the stress: [ts, tn]
                this.intPoint(i).stressOld = [stressNM(3);stressNM(2)];
                this.intPoint(i).stress    = [stressNM(3);stressNM(2)];

            end

        end

        function cohesiveStress = geostaticStressAtPoint(this,X, elem, K0, yTop)

            % Specific weight
                rhob = elem.intPoint(1).constitutiveMdl.rhoBulk();
                grav = elem.intPoint(1).constitutiveMdl.gravityAcc();
                gw   = elem.intPoint(1).constitutiveMdl.waterSpecificWeight();

                % Compute the effective stress components in the global
                % cartesian coordinate system
                sy = -100-(yTop - X(2))*(grav*rhob - gw);
                sx = K0*sy;

                stressXY = [sx;sy;0.0];

                stressNM = this.rotXY2nm*stressXY;

                cohesiveStress = [stressNM(3);stressNM(2)];

        end

        %------------------------------------------------------------------
        % Returns the fracture element matrices and vectors
        function [Kf, Lc, Ldc, Ld, Hf, Qf, Sf, Ef, fi, Ks] = elementData(this,Ue,pf,idFract,elem)

            % Initialize the sub-matrices of the fracture
            Kf      = zeros(this.nglu,this.nglu);  
            Ks      = zeros(this.nglu,elem.nglutot);
            Hf      = zeros(this.ndofPf,this.ndofPf);
            Qf      = zeros(elem.ngluenr,this.ndofPf);
            Sf      = zeros(this.ndofPf,this.ndofPf);
            Ef      = zeros(this.ndofPf,elem.nglutot);
            Lc_pp   = zeros(elem.nglp,elem.nglp);
            Lc_pDp  = zeros(elem.nglp,elem.nglpenrDp);
            Lc_DpDp = zeros(elem.nglpenrDp,elem.nglpenrDp);
            Ldc_ppf = zeros(elem.nglp,this.ndofPf);
            Ldc_apf = zeros(elem.nglpenrDp,this.ndofPf);
            Ld      = zeros(this.ndofPf,this.ndofPf);
             
            % Initialize the mechanical part of the internal force vector
            fi = zeros(this.nglu, 1);
            
            % Compute the rotation matrix: [x y] => [shear normal]
            R = this.rotationPointMtrx();

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints

                % Numerical integration term. The determinant is ld/2.
                c = this.intPoint(i).w * this.ld/2 * this.t;

                % Cartesian coordinates of the integration point 
                X = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Natural coordinates associated with the continuum element
                % of this point
                Xn = elem.shape.coordCartesianToNatural(elem.node,X);

                % Shape function matrix of the continuum
                N = elem.shape.shapeFncMtrx(Xn);

                % --- Mechanical part -------------------------------------

                % Shape function matrix
                Nu = this.interpJumpShapeMtrx(this.intPoint(i).X,elem.enrVar);

                % Evaluate the jump at the integration point in the local
                % coordinate system
                w = R * Nu * Ue;
           
                % Compute the stress vector and the constitutive matrix
                [td,T] = this.intPoint(i).constitutiveModel(w);

                % Shape function matrix of the fracture to interpolate pf
                Nf = this.shapeFncMtrx(this.intPoint(i).X);

                % Discontinuity pressure at the integration point
                pfIP = Nf * pf;

                % Get Biot's coefficient
                biot = this.intPoint(i).constitutiveMdl.biotCoeff();

                % Numerical integration of the stiffness matrix and the
                % internal force vector
                Kf = Kf + Nu' * R' *  T * R * Nu * c;
                fi = fi + Nu' * (R' * td - this.n' * biot * pfIP) * c;

                % --- Hydraulic part --------------------------------------
                
                % Enhanced shape function matrix
                Ntop = elem.topEnhancedShapeFncMtrx(N,Xn,idFract);
                Nbot = elem.bottomEnhancedShapeFncMtrx(N,Xn,idFract);

                % Gradient of the shape function matrix
                dNfds = this.gradShapeFncMtrx(this.intPoint(i).X);
               
                % Compute longitudinal permeability
                kl = this.intPoint(i).constitutiveMdl.computeLongPermeability(this.intPoint(i));
        
                % Get transverse flow leak-off parameters
                [cT, cB] = this.intPoint(i).constitutiveMdl.getLeakOffParameters(this.intPoint(i));

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibilityCoeff(this.intPoint(i));

                % Numerical integration of the coupling matrix associated
                % with the porous-media dofs
                Lc_pp   = Lc_pp   + (N' * N) * (cT + cB) * c;
                Lc_pDp  = Lc_pDp  + ((N' * Ntop) * cT + (N' * Nbot) * cB) * c;
                Lc_DpDp = Lc_DpDp + ((Ntop' * Ntop) * cT + (Nbot' * Nbot) * cB) * c;

                % Numerical integration of the coupling matrix associated
                % with the discontinuity dofs
                Ldc_ppf = Ldc_ppf + (N' * Nf) * (cT + cB)  * c;
                Ldc_apf = Ldc_apf + ((Ntop' * Nf) * cT + (Nbot' * Nf) * cB) * c;
        
                % Numerical integration of the fracture fluid-flow matrices
                Hf = Hf  + (dNfds' * dNfds) * kl * c;
                Ld = Ld  + (Nf' * Nf) * (cT + cB)  * c;

                % Numerical integration of the fracture compressibility
                % matrix
                Sf = Sf + (Nf' * Nf) * comp * c;

                % --- HM-Coupling part ------------------------------------
                
                % Numerical integration of the coupling matrix
                Qf = Qf + Nu' * this.n' * biot * Nf * c;
                % 
                % % Enrichment mapping matrix for the displacement field
                % Mu = elem.elementDisplacementMappingMtrx();
                % 
                % % Compute the B matrix at the int. point and the detJ
                % [dNdx, ~] = elem.shape.dNdxMatrix(elem.node,Xn);
                % 
                % % Assemble the B-matrix for the mechanical part
                % Bu = elem.shape.BMatrix(dNdx);
                % 
                % % Compute the enriched B-matrix
                % BuEnrMid = elem.enhancedMiddleCompatibilityMtrx(dNdx,Xn);
                % 
                % % Assemble the B-matrix for the mechanical part
                % BuEnrMid = elem.shape.BMatrix(BuEnrMid);
                % BuEnrMid = BuEnrMid * Mu;
                % Bmidaug  = [Bu , BuEnrMid];
                % 
                % wn = this.intPoint(i).constitutiveMdl.getAperture(this.intPoint(i));
                % Ef = Ef + Nf' * biot * wn * this.rs * Bmidaug  * c;

            end

            % Assemble the flow coupling matrices
            Lc = [ Lc_pp ,  Lc_pDp;
                   Lc_pDp', Lc_DpDp];

            Ldc = [ Ldc_ppf; Ldc_apf ];
            
        end

        % -----------------------------------------------------------------
        % Update the traction stress vector given a displacement jump in
        % the global cartesian system
        function updateTractionStress(this,w)
           
            % Compute the stress vector and the constitutive matrix
            this.intPoint(1).constitutiveModel(w);

        end

        % -----------------------------------------------------------------
        % Get fracture condensation matrix of the longitudinal pressure dof
        function T = condensationMtrx(this,elem,idfrac)

            % Natural coordinates associated with the continuum element
            % of nodes of the fracture
            Xn1 = elem.shape.coordCartesianToNatural(elem.node,this.node(1,:));
            Xn2 = elem.shape.coordCartesianToNatural(elem.node,this.node(2,:));

            % Shape function matrix of the continuum
            N1 = elem.shape.shapeFncMtrx(Xn1);
            N2 = elem.shape.shapeFncMtrx(Xn2);

            Nenr1 = elem.middleEnhancedShapeFncMtrx(N1,Xn1,idfrac);
            Nenr2 = elem.middleEnhancedShapeFncMtrx(N2,Xn2,idfrac);

            % Assemble the condensation matrix
            T = [N1 Nenr1 ; N2 Nenr2];
        end

        %------------------------------------------------------------------
        function [td, Td, S, DBaug] = stabilizedCohesiveLaw(this, idFract, ip, elem, td, Td)

            % Get stabilization matrix
            S = this.stabilizationMtrx(ip,elem);

            % Cartesian coordinates of the integration point 
            X = this.shape.coordNaturalToCartesian(this.node,ip.X);

            % Compute the weighted Cauchy stress vector
            [stress_gamma, DBaug] = elem.weightedStress(idFract,X);

            % Rotate stress vector to the discontinuity local system
            stress_gamma = this.rotXY2nm*stress_gamma;

            % Traction vector [ts; tn]
            t_gamma = [stress_gamma(3); stress_gamma(2)];

            % Stabilized cohesive stress vector
            td = S * td - (eye(2) - S) * t_gamma;

            % Stabilized constitutive matrix
            Td = S * Td;
        end

        %------------------------------------------------------------------
        % This function computes the tangential local coordinate s for a
        % given coordinate in the cartesian system.
        function s = tangentialLocCoordinate(this,X)

            % Relative position vector
            DX = X - this.Xref;

            % Tangential coordinate
            s = this.m*DX';
            if abs(s) > this.ld/2.0
                s = sign(s)*this.ld/2.0;
            end

        end

        %------------------------------------------------------------------
        % Compute the cartesian coordinates of a point given the tangential
        % coordinate s.
        function X = cartesianCoordinate(this,s)

            % Cartesian coordinates of the given point
            X = this.Xref + s*this.m;

        end

        %------------------------------------------------------------------
        % Compute the stabilization matrix
        function S = stabilizationMtrx(this,ip,elem)

            % Initialize the stabilization matrix
            S = eye(2,2);

            % Get the stabilization parameters
            % beta = this.stabilizationParameter(elem);
            % 
            % % Get the initial stiffness parameters
            % ks = ip.constitutiveMdl.getShearStiffness();
            % kn = ip.constitutiveMdl.getNormalStiffness();
            % 
            % % Compute the shear stabilization component
            % S(1,1) = beta / (ks + beta);
            % 
            % % Compute the normal stabilization component
            % S(2,2) = beta / (kn + beta);

        end

        %------------------------------------------------------------------
        % Compute the stabilization parameter
        function beta = stabilizationParameter(this,elem)

            % Get the area of the element
            A = elem.getDomainArea();

            % Get the area of the sub-region formed by this fracture
            Atop = elem.areaPositiveSubRegion(this);

            % Get the Young modulus of the continuum
            E = elem.intPoint(1).constitutiveMdl.youngModulus();

            % Compute the stabilization parameter
            beta = 2.0 * E * (Atop/(A*A) + (A - Atop)/(A*A)) * this.ld;

        end

        %------------------------------------------------------------------
        % This function compute the stress interpolation vector
        function S = stressIntVctFnc(this, shape, node, intpOrder)

            % Initialize the Gram matrix
            S = zeros(shape.getSizeStressIntVct(), intpOrder+1);
 
            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints

                % Tangential coordinate and point in the global coordinate
                % system
                s = this.intPoint(i).X * this.ld/2.0;

                % Relative position coordinate
                Xrel = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Compute the integrand of the stress interpolation vector
                dS = shape.integrandStressIntVct(s,Xrel,intpOrder);
        
                % Numerical integration term. The determinant is ld/2.
                c = this.intPoint(i).w * this.ld/2;
        
                % Numerical integration of the stiffness matrix and the
                % internal force vector
                S = S + dS * c;

            end
            % S = S / this.ld;
        end

        function P = stressProjectionMatrix(this)

            P = [this.n(1) , 0.0 ;
                 0.0       , this.n(2);
                 this.n(2) , this.n(1)];

        end

    end

end