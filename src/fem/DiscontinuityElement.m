%% DiscontinuityElement Class
% This in an abstract class that defines a discontinuity element in a finite element mesh.
% It provides methods to compute geometric and physical properties of the discontinuity.
% 
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% 
%% Class definition
classdef DiscontinuityElement < handle    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        shape      = [];    % Object of the Shape class
        node       = [];    % Nodes of the fem mesh
        t          = 1.0;   % Thickness
        mat        = [];    % Material object
        intOrder   = 2;     % Order of the numerical integration
        dof        = [];    % Degrees of freedom vector
        dofOld     = [];    % Old degrees of freedom
        ndof       = 0;     % Number of dofs
        nIntPoints = 1;     % Number of integration points
        intPoint   = [];    % Vector with integration point objects       
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement(node,mat)
            if (nargin > 0)
                this.shape = Shape_Bar();
                this.node  = node;
                this.mat   = mat;
            end
        end
    end

    %% Abstract methods
    methods(Abstract)
        %------------------------------------------------------------------
        % Assemble element matrices and vectors.
        % Outputs:
        %    Ke : element "stiffness" matrix
        %    Ce : element "damping" matrix
        %    fe : element "external force" vector
        %    fi : element "internal force" vector
        % dfidu : element matrix of derivative of the internal force wrt displacement
        [Ke,Ce,fi,fe,dfidu] = elementData(this,ae);
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Compute discontinuity reference point.
        function Xr = referencePoint(this)
            Xr = 0.5 * (this.node(1,:) + this.node(2,:));
        end

        %------------------------------------------------------------------
        % Compute discontinuity length.
        function l = ld(this)
            dx = this.node(2,1) - this.node(1,1);
            dy = this.node(2,2) - this.node(1,2);
            l  = sqrt(dx.^2 + dy.^2);
        end

        %------------------------------------------------------------------
        % Compute discontinuity tangential vector.
        function m = tangentialVector(this)
            DX = this.node(2,:) - this.node(1,:);
            m  = DX' / norm(DX);
        end

        %------------------------------------------------------------------
        % Compute discontinuity normal vector.
        % Defined considering n = ez * m, where ez = [0 0 1].
        function n = normalVector(this)
            m = this.tangentialVector();
            n = [-m(2) ; m(1)];
        end

        %------------------------------------------------------------------
        % Compute heaviside function associated with the discontinuity at a given point.
        function h = heaviside(this,X)
            n  = this.normalVector();
            Xr = this.referencePoint();
            DX = X - Xr;
            h  = max(sign(DX*n),0.0);
        end

        %------------------------------------------------------------------
        % Initialize degrees of freedom vector.
        function initializeDofs(this,ndofs)
            this.dof = 1:this.ndof;
            this.dof = this.dof + ndofs;
        end

        %------------------------------------------------------------------
        % Update state variables.
        function updateStateVar(this)
            this.dofOld = this.dof;
            for i = 1:this.nIntPoints
                this.intPoint(i).updateStateVar();
                this.intPoint(i).updateStressVct();
                this.intPoint(i).updateStrainVct();
            end
        end
    end
end
