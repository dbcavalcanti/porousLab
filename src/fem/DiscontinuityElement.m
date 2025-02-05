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
        t            = 1.0;           % Thickness
        mat          = [];            % Material object
        intOrder     = 2;             % Order of the numerical integration
        ndof         = 1;             % Number of dofs
        nIntPoints   = 1;             % Number of integration points
        intPoint     = [];            % Vector with integration point objects       
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement(node, mat)
            if (nargin > 0)
                this.shape   = Shape_Bar();
                this.node    = node;
                this.mat     = mat;
            end
        end
    end
    %% Abstract methods
    methods(Abstract)
        [Ke, Ce, fi, fe, dfidu] = elementData(this, ae);
    end
    %% Public methods
    methods
        % -----------------------------------------------------------------
        % Compute the discontinuity reference point
        function Xr = referencePoint(this)
            Xr = 0.5*(this.node(1,:) + this.node(2,:));
        end

        % -----------------------------------------------------------------
        % Compute the discontinuity length
        function l = ld(this)
            dx = this.node(2,1) - this.node(1,1);
            dy = this.node(2,2) - this.node(1,2);
            l  = sqrt(dx.^2 + dy.^2);
        end

        % -----------------------------------------------------------------
        % Compute the discontinuity tangential vector
        function m = tangentialVector(this)
            DX = this.node(2,:) - this.node(1,:);
            m  = DX'/norm(DX);
        end

        % -----------------------------------------------------------------
        % Compute the discontinuity normal vector
        % Defined considering n = ez x m, where ez = [0 0 1]
        function n = normalVector(this)
            m = this.tangentialVector();
            n = [-m(2) ; m(1)];
        end

        % -----------------------------------------------------------------
        % Compute the heaviside function associated with this discontinuity
        function h = heaviside(this,X)
            n  = this.normalVector();
            Xr = this.referencePoint();
            DX = X - Xr;
            h  = max(sign(DX*n),0.0);
        end

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