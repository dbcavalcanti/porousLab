%% Element class
%
% This class defines a finite element model (consider a ISOQ4 element)
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: December 2022
%%%
% Initially prepared for the course CIV 2801 - Fundamentos de Computação
% Gráfica, 2022, second term, Department of Civil Engineering, PUC-Rio.
%
%% Class definition
classdef RegularElement < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        type       = 'ISOQ4';       % type of element
        shape      = [];            % Object of the Shape class
        node       = [];            % Nodes of the fem mesh
        connect    = [];            % Nodes connectivity
        t          = 1.0;           % Thickness
        mat        = [];            % Vector with material properties
        intOrder   = 2;             % Order of the numerical integration
        nnd_el     = 4;             % Number of nodes per element
        ndof_nd    = 1;             % Number of dof per node
        glp        = [];            
        glpg       = [];            % Vector of the regular degrees of freedom
        gle        = [];            % Vector of the degrees of freedom
        nglp       = 0;             % Number of regular p-dof
        ngle       = 0;             % Number of total dof
        ue         = [];            % Element's displacement vector
        due        = [];            % Element's increment displacement
        nIntPoints = 1;             % Number of integration points
        intPoint   = [];            % Vector with integration point objects
        result     = [];            % Result object to plot the results
        isEnriched = false;         % Flag to check if the element is enriched
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement(type, node, elem, t, ...
                mat, intOrder, glp, glpg)
            if (nargin > 0)
                if strcmp(type,'ISOQ4')
                    this.shape = Shape_ISOQ4();
                elseif strcmp(type,'ISOQ8')
                    this.shape = Shape_ISOQ8();
                elseif strcmp(type,'CST')
                    this.shape = Shape_CST();   
                elseif strcmp(type,'LST')
                    this.shape = Shape_LST();
                end
                this.node     = node;
                this.nnd_el   = size(node,1);
                this.connect  = elem;
                this.t        = t;
                this.mat      = mat;
                this.intOrder = intOrder;
                this.glp      = glp;
                this.glpg     = glpg;
                this.gle      = [glp , glpg];
                if (length(this.glp) ~= length(this.glpg))
                    error('Wrong number of pressure dofs');
                end
                this.nglp     = length(this.glp);
                this.ngle     = length(this.gle);
                order = this.sortCounterClockWise(this.node);
                this.result   = Result(this.node(order,:),1:length(this.connect),0.0*ones(this.nnd_el,1),'Model');
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(this.intOrder);

            % Initialize the integration points objects
            
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                constModel = MaterialTwoPhaseFlow(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
            end
            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        % This function assembles the element matrices and vectors 
        %
        % Output:
        %   Ke : element "stiffness" matrix
        %   Ce : element "damping" matrix
        %   fe : element "internal force" vector
        %
        function [Ke, Ce, fe] = elementData(this)

            % Initialize the sub-matrices
            Hll = zeros(this.nglp, this.nglp);
            Hgg = zeros(this.nglp, this.nglp);  
            Sll = zeros(this.nglp, this.nglp);
            Sgg = zeros(this.nglp, this.nglp);
            Sgl = zeros(this.nglp, this.nglp);
            Slg = zeros(this.nglp, this.nglp);
            O   = zeros(this.nglp, this.nglp);
            
            % Vector of the nodal pore-pressure dofs
            pl = this.ue(1:this.nglp);
            pg = this.ue(1+this.nglp:end);

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Capillary pressure at the integration point
                pcIP = Np * (pg - pl);

                % Compute the saturation degree at the integration point
                Sl = this.intPoint(i).constitutiveMdl.saturationDegree(pcIP);
        
                % Compute the permeability matrix
                [kl, kg] = this.intPoint(i).constitutiveMdl.permeabilityMtrcs(Sl);

                % Get compressibility coefficients
                [cll,cgg,clg,cgl] = this.intPoint(i).constitutiveMdl.compressibilityCoeffs(pcIP,Sl);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
        
                % Compute permeability sub-matrices
                Hll = Hll + Bp' * kl * Bp * c;
                Hgg = Hgg + Bp' * kg * Bp * c;

                % Compute compressibility matrices
                Sll = Sll + Np' * cll * Np * c;
                Sgg = Sgg + Np' * cgg * Np * c;
                Slg = Slg + Np' * clg * Np * c;
                Sgl = Sgl + Np' * cgl * Np * c;
            end

            % Assemble the element matrices
            Ke = [ Hll, O;
                    O , Hgg ];

            Ce = [ Sll , Slg;
                   Sgl , Sgg ];


            % Assemble element internal force vector
            fe = zeros(this.ngle, 1);
            fe(1:this.nglp)     = Hll*pl;
            fe(1+this.nglp:end) = Hgg*pg;
            
        end

        %------------------------------------------------------------------
        % Returns the internal force vector
        function fe = elementInternalForceVector(this,Ue)

            fe = zeros(this.ngle, 1);
            fi = zeros(this.nglp, 1);

            % Initialize Fluid-flow matrix
            H = zeros(this.nglpg, this.nglpg);    

            % Initialize 2D identity vector
            m = [1.0 ; 1.0 ; 0.0];

            % Vector of the nodal pore-pressure dofs
            pe = Ue(1+this.nglp:end);

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints
                
                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Pore pressure at the integration point
                pIP = Np * pe;

                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.shape.BMatrix(Bp);

                % Compute the increment of the strain vector
                strain = Bu*Ue(1:this.nglp);
        
                % Compute the stress vector and the constitutive matrix
                [stress,~] = this.intPoint(i).constitutiveModel(strain); 

                % Compute the permeability matrix
                Dh = this.intPoint(i).constitutiveMdl.permeabilityMtrx();

                % Get Biot's coefficient
                biot = this.intPoint(i).constitutiveMdl.biotCoeff();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t; 

                % Fluid-flow matrix
                H = H + Bp' *  Dh  * Bp * c;

                % Numerical integration of the internal force vector
                fi = fi + Bu' * (stress - biot * pIP * m) * c;

            end

            % Hydraulic part of the internal force vector
            fp = H*pe;

            % Assemble element internal force vector
            fe(1:this.nglp)     = fi;
            fe(1+this.nglp:end) = fp;
            
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

            % Get nodal pressures
            pl = this.ue(1:this.nglp);

            % capillary field
            p = Nm*pl;
        
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = gasPressureField(this,X,ue)
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

            % Get nodal pressures
            pg = this.ue(1+this.nglp:end);

            % capillary field
            p = Nm*pg;
        
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = capillaryPressureField(this,X,ue)
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

            % Get nodal pressures
            pl = this.ue(1:this.nglp);
            pg = this.ue(1+this.nglp:end);

            % capillary field
            p = Nm*(pg - pl);
        
        end

        %------------------------------------------------------------------
        % Function to compute the body force vector of the element.
        function Fb = bodyForce(this)
        
            Fb = zeros(this.nglp, 1);

            grav = [0.0; -1.0];

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                N = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Compute the B matrix at the int. point and the detJ
                [~, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                Nu = this.shape.NuMtrx(N);
                
                % Specific mass of the bulk
                rhoBulk = this.intPoint(i).constitutiveMdl.rhoBulk();
                gAcc = this.intPoint(i).constitutiveMdl.gravityAcc();

                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;

                % Numerical integration of the internal force vector
                Fb = Fb + Nu' *grav * rhoBulk * gAcc * c;

            end
        
        end

        % -----------------------------------------------------------------
        % Function to update the state variables
        function updateStateVar(this)

            for i = 1:this.nIntPoints
                this.intPoint(i).updateStateVar();
            end

        end

        %------------------------------------------------------------------
        % Function to compute the area of the element domain
        function A = getDomainArea(this)
            A = this.calculateArea(this.node);
        end

        function A = calculateArea(~,node)

            % Vertices of the coordinates
            vx = node(:,1); 
            vy = node(:,2);
        
            % Shifted vertices
            vxS = vx([2:end 1]);
            vyS = vy([2:end 1]); 
        
            % Compute the area of the polygon
            temp = vx.*vyS - vy.*vxS;
            A    = 0.5*sum(temp);
            
        end
        

    end

    %% Public methods associated with the pos-processing
    methods

        %------------------------------------------------------------------
        % Update the result's object vertices property
        % If the 'Undeformed' configuration is selected, nothing needs to
        % be done.
        function updateResultVertices(this,configuration)
            if strcmp(configuration,'Deformed')
                Nodes = this.getDeformedConfiguration();
                this.result.setVertices(Nodes);
            end  
        end

        %------------------------------------------------------------------
        % Update the result's object vertices property
        % If the 'Undeformed' configuration is selected, nothing needs to
        % be done.
        function updateResultFaces(this,faces)
            this.result.setFaces(faces);
        end
        
        %------------------------------------------------------------------
        % Update the result's object vertices data property
        function updateResultVertexData(this,type)
            this.result.setDataLabel(type);
            switch type
                case 'Ux'
                    ndResults = this.getNodalDisplacementField(1);
                case 'Uy'
                    ndResults = this.getNodalDisplacementField(2);
                case 'Sxx'
                    ndResults = this.getNodalStressField(1);
                case 'Syy'
                    ndResults = this.getNodalStressField(2);
                case 'Sxy'
                    ndResults = this.getNodalStressField(3);
            end
            this.result.setVertexData(ndResults);
        end

    end
    methods(Static)

        %------------------------------------------------------------------
        % This function sorts counterclockwise a set of nodes.
        % It uses as a reference point the centroid defined by the nodes.
        function order = sortCounterClockWise(NODE)
            
            % Centroid coordinate
            cx = mean(NODE(:,1));
            cy = mean(NODE(:,2));
            
            % Compute the angle that the relative vector of the vertices 
            % from the centroid has with the horizontal axis
            a = atan2(NODE(:,2) - cy, NODE(:,1) - cx);
            
            % Sort the angles
            [~, order] = sort(a);
            
        end

    end
end