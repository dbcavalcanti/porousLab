%% RegularElement class
%
% This class defines a finite element of a two-phase flow formulation using
% the liquid pressure (Pl) and the gas pressure (Pg) as primary
% variables.
%
%% Author
% Danilo Cavalcanti
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
        massLumping = false;        % Flag to apply a diagonalization of the compressibility matrix
        lumpStrategy = 1;
        isAxisSymmetric = false;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement(type, node, elem, t, ...
                mat, intOrder, glp, glpg, massLumping, lumpStrategy,isAxisSymmetric)
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
                this.massLumping = massLumping;
                this.lumpStrategy = lumpStrategy;
                this.isAxisSymmetric = isAxisSymmetric;
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
        function [Ke, Ce, fi, fe] = elementData(this)

            % Initialize the sub-matrices
            Hll = zeros(this.nglp, this.nglp);
            Hlg = zeros(this.nglp, this.nglp);
            Hgl = zeros(this.nglp, this.nglp);
            Hgg = zeros(this.nglp, this.nglp);  
            Sll = zeros(this.nglp, this.nglp);
            Slg = zeros(this.nglp, this.nglp);
            Sgl = zeros(this.nglp, this.nglp);
            Sgg = zeros(this.nglp, this.nglp);

            % Initialize external force vector
            fel = zeros(this.nglp, 1);
            feg = zeros(this.nglp, 1);
            
            % Vector of the nodal pore-pressure dofs
            pc = this.getNodalCapillaryPressure();
            pg = this.getNodalGasPressure();

            % Initialize the volume of the element
            vol = 0.0;

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Pressure values at the integration point
                pcIP = Np * pc;
                pgIP = Np * pg;

                % Compute the saturation degree at the integration point
                Sl = this.intPoint(i).constitutiveMdl.saturationDegree(pcIP);
        
                % Compute the permeability matrix
                [kll, klg, kgl, kgg] = this.permeabilityTensors(this.intPoint(i),pgIP,pcIP,Sl);

                % Get compressibility coefficients
                [cll, clg, cgl, cgg] = this.compressibilityCoeffs(this.intPoint(i),pgIP,pcIP,Sl);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
        
                % Compute permeability sub-matrices
                Hll = Hll + Bp' * kll * Bp * c;
                Hlg = Hlg + Bp' * klg * Bp * c;
                Hgl = Hgl + Bp' * kgl * Bp * c;
                Hgg = Hgg + Bp' * kgg * Bp * c;

                % Compute compressibility matrices
                if ((this.massLumping) && (this.lumpStrategy == 1))
                    Sll = Sll + diag(cll*Np*c);
                    Sgg = Sgg + diag(cgg*Np*c);
                    Slg = Slg + diag(clg*Np*c);
                    Sgl = Sgl + diag(cgl*Np*c);
                elseif (this.massLumping == false)
                    Sll = Sll + Np' * cll * Np * c;
                    Sgg = Sgg + Np' * cgg * Np * c;
                    Slg = Slg + Np' * clg * Np * c;
                    Sgl = Sgl + Np' * cgl * Np * c;
                end
                
                % Compute the gravity forces
                if (this.mat.porousMedia.gravityOn)
                    [fel,feg] = this.addGravityForces(fel,feg,Bp,kll,kgg,pgIP-pcIP,pgIP,c);
                end

                % Compute the element volume
                vol = vol + c;
            end

            % Compute the lumped mass matrix
            if ((this.massLumping) && (this.lumpStrategy == 2))
                [Sll,Slg,Sgl,Sgg] = lumpedCompressibilityMatrices(this, pc, pg, vol);
            end

            % Assemble the element permeability
            Ke = [ Hll, Hlg;
                   Hgl, Hgg ];

            % Assemble the element compressibility matrix
            Ce = [ Sll , Slg;
                   Sgl , Sgg ];

            % Assemble element internal force vector
            fi = Ke * this.ue;

            % Assemble element external force vector
            fe = [fel; feg];
            
        end

        % -----------------------------------------------------------------
        % Compute the permeability tensors
        function [kll, klg, kgl, kgg] = permeabilityTensors(~,ip,pg,pc,Sl)
             [kll, klg, kgl, kgg] = ip.constitutiveMdl.permeabilityMtrcs(Sl,pg-pc,pg);
        end

        % -----------------------------------------------------------------
        % Compute the compressibility coefficients
        function [cll, clg, cgl, cgg] = compressibilityCoeffs(~,ip,pg,pc,Sl)
             [cll, clg, cgl, cgg] =  ip.constitutiveMdl.compressibilityCoeffs(Sl,pg-pc,pg);
        end

        %------------------------------------------------------------------
        % Compute the lumped mass matrices
        function [Sll,Slg,Sgl,Sgg] = lumpedCompressibilityMatrices(this, pc, pg, vol)

            % Shape function matrix
            Np = this.shape.shapeFncMtrx([0.0,0.0]);

            % Pressure values at the integration point
            pcIP = Np * pc;
            pgIP = Np * pg;

            % Compute the saturation degree at the integration point
            Sl = this.intPoint(1).constitutiveMdl.saturationDegree(pcIP);

            % Get compressibility coefficients
            [cll, clg, cgl, cgg] = this.compressibilityCoeffs(this.intPoint(1),pgIP,pcIP,Sl);

            % Mass distribution factor
            factor = vol / this.nnd_el;

            % Compressibility matrices
            Sll = cll * factor * eye(this.nglp,this.nglp);
            Slg = clg * factor * eye(this.nglp,this.nglp);
            Sgl = cgl * factor * eye(this.nglp,this.nglp);
            Sgg = cgg * factor * eye(this.nglp,this.nglp);

        end

        %------------------------------------------------------------------
        % Add contribution of the gravity forces to the external force vct
        function [fel,feg] = addGravityForces(this,fel,feg,Bp,kl,kg,pl,pg,c)

            % Get gravity vector
            grav = this.mat.porousMedia.g * this.mat.porousMedia.b;

            % Get fluid densities
            rhol = this.mat.liquidFluid.getDensity();
            rhog = this.mat.gasFluid.getDensity();

            % Compute the contribution of the gravitational forces
            fel = fel + Bp' * kl * rhol * grav * c;
            feg = feg + Bp' * kg * rhog * grav * c;
            
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the liquid pressure
        function pl = getNodalLiquidPressure(this)
            pl = this.ue(1:this.nglp);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the gas pressure
        function pg = getNodalGasPressure(this)
            pg = this.ue(1+this.nglp:end);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the capillary pressure
        function pc = getNodalCapillaryPressure(this)
            pl = this.getNodalLiquidPressure();
            pg = this.getNodalGasPressure();
            pc = pg - pl;
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
            pl = this.getNodalLiquidPressure();

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
            pg = this.getNodalGasPressure();

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
            pc = this.getNodalCapillaryPressure();

            % capillary field
            p = Nm*pc;
        
        end

        %------------------------------------------------------------------
        % Function to compute the liquid saturation field inside a given element
        function Sl = liquidSaturationField(this,X,ue)

            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Capillary pressure at the given point
            pc = Nm*this.getNodalCapillaryPressure();

            % Compute the liquid saturation degree
            Sl = this.intPoint(1).constitutiveMdl.saturationDegree(pc);
        
        end

        %------------------------------------------------------------------
        % Function to compute the gas saturation field inside a given element
        function Sg = gasSaturationField(this,X,ue)
            
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Capillary pressure at the given point
            pc = Nm*this.getNodalCapillaryPressure();

            % Compute the liquid saturation degree
            Sl = this.intPoint(1).constitutiveMdl.saturationDegree(pc);

            % Gas saturation degree
            Sg = 1.0 - Sl;
        
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