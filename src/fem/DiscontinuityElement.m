%% RegularElement class
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef DiscontinuityElement < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
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
        ueOld        = [];            % Element's old displacement vector
        due          = [];            % Element's increment displacement
        nIntPoints   = 1;             % Number of integration points
        intPoint     = [];            % Vector with integration point objects
        massLumping  = false;         % Flag to apply a diagonalization of the compressibility matrix
        lumpStrategy = 1;             % Id of the diagonalization strategy
        DTime        = [];            % Time increment          
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement(node, elem, t, ...
                mat, intOrder, massLumping, lumpStrategy)
            if (nargin > 0)
                this.shape = Shape_Bar();
                this.node     = node;
                this.nnd_el   = size(node,1);
                this.connect  = elem;
                this.t        = t;
                this.mat      = mat;
                this.intOrder = intOrder;
                this.massLumping = massLumping;
                this.lumpStrategy = lumpStrategy;
            end
        end
    end

    %% Abstract methods
    methods(Abstract)
        [Ke, Ce, fi, fe, dfidu] = elementData(this);
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
    end
end