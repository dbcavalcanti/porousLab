%% Mechanical discontinuity element class
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef DiscontinuityElement_M < DiscontinuityElement    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        stretchingMode  = false;
        relRotationMode = false;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement_M(node, mat)
            this = this@DiscontinuityElement(node, mat)
            this.ndof = 2;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        function addStretchingMode(this,flag)
            this.stretchingMode = flag;
            if flag == true
                this.ndof = this.ndof + 1;
            end
        end

        %------------------------------------------------------------------
        function addRelRotationMode(this,flag)
            this.relRotationMode = flag;
            if flag == true
                this.ndof = this.ndof + 1;
            end
        end

        %------------------------------------------------------------------
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(1);

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                constModel = MaterialDiscontinuity_M(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
                intPts(i).initializeMechanicalAnalysisModel('Interface');
            end
            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        function [Ke, Ce, fi, fe, dfidu] = elementData(this, ae)
            
            % Declare output matrices that won't be used
            Ce = []; fe = []; dfidu = [];

            % Initialize the matrices for the numerical integration
            Ke = zeros(this.ndof,this.ndof);
            fi = zeros(this.ndof,1);

            % Get the lenght of the discontinuity
            ld = this.ld();

            % Get the discontinuity reference point
            Xr = this.referencePoint();

            % Get the discontinuity tangential vector
            m = this.tangentialVector();

            % Initialize output matrices
            for i = 1:this.nIntPoints

                % Cartesian coordinates of the integration point 
                X = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Get the shape function matrix
                Nd = this.enrichmentInterpolationMatrix(X,Xr,m);

                % Compute the strain vector
                this.intPoint(i).strain = Nd * ae;

                % Compute the stress vector and the constitutive matrix
                [td,Td] = this.intPoint(i).mechanicalLaw();

                % Numerical integration term. The determinant is ld/2.
                c = 0.5 * ld * this.intPoint(i).w * this.t;

                % Compute the stiffness sub-matrix
                Ke = Ke + Nd' * Td * Nd * c;

                % Compute the internal force vector
                fi = fi + Nd' * td * c;

            end
        end

        %------------------------------------------------------------------
        function Nd = enrichmentInterpolationMatrix(this,X,Xr,m)
            Nd = zeros(2,this.ndof);
            Nd(1,1) = 1.0;
            Nd(2,2) = 1.0;
            if this.ndof > 2
                s = m' * (X' - Xr');
                c = 3;
                if this.stretchingMode
                    Nd(1,c) = s;
                    c = c + 1;
                end
                if this.relRotationMode
                    Nd(2,c) = s;
                end
            end
        end
    end
end