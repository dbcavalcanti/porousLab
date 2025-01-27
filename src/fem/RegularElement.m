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
        type         = 'ISOQ4';       % type of element
        shape        = [];            % Object of the Shape class
        node         = [];            % Nodes of the fem mesh
        connect      = [];            % Nodes connectivity
        t            = 1.0;           % Thickness
        mat          = [];            % Vector with material properties
        intOrder     = 2;             % Order of the numerical integration
        nnd_el       = 4;             % Number of nodes per element
        ndof_nd      = 1;             % Number of dof per node
        gle          = [];            % Vector of the degrees of freedom
        ngle         = 0;             % Total number of dofs
        ue           = [];            % Element's displacement vector
        due          = [];            % Element's increment displacement
        nIntPoints   = 1;             % Number of integration points
        intPoint     = [];            % Vector with integration point objects
        result       = [];            % Result object to plot the results
        isEnriched   = false;         % Flag to check if the element is enriched
        massLumping  = false;         % Flag to apply a diagonalization of the compressibility matrix
        lumpStrategy = 1;             % Id of the diagonalization strategy
        isAxisSymmetric = false;      % Flag to axissymetric models
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement(type, node, elem, t, ...
                mat, intOrder, massLumping, lumpStrategy,isAxisSymmetric)
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

        % -----------------------------------------------------------------
        % Function to update the state variables
        function updateStateVar(this)

            for i = 1:this.nIntPoints
                this.intPoint(i).updateStateVar();
                this.intPoint(i).updateStressVct();
                this.intPoint(i).updateStrainVct();
            end

        end

        %------------------------------------------------------------------
        % Function to compute the element characteristic length
        function lc = characteristicLength(this)
            lc = this.getDomainArea();
            if strcmp(this.shape.type,'CST') || strcmp(this.shape.type,'LST')
                lc = lc * sqrt(2.0);
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