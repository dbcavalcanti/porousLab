%% DiscontinuityElement_M Class
% This class defines a mechanical discontinuity element that extends the
% _DiscontinuityElement_ base class. It includes additional functionality
% for handling stretching and relative rotation modes, as well as methods
% for initializing integration points and computing element data.
%
%% Methods
% * *addStretchingMode*: Enables or disables the stretching mode based on 
%                        the flag.
% * *addRelRotationMode*: Enables or disables the relative rotation mode 
%                         based on the flag.
% * *initializeIntPoints*: Initializes the integration points for the 
%                          discontinuity element. Retrieves integration 
%                          points' coordinates and weights, and creates 
%                          _IntPoint_ objects with the associated material 
%                          model.
% * *elementData*: Computes the element stiffness matrix, internal force 
%                  vector, and other element data based on the input 
%                  displacement vector. Performs numerical integration 
%                  over the integration points.
% * *enrichmentInterpolationMatrix*: Computes the enrichment interpolation 
%                                    matrix for the given integration point
%                                    coordinates, reference point, and 
%                                    tangential vector.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00: Initial version (January 2024).
%
%% Class Definition
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
        % Enables the stretching mode. If enables, the number of degrees of
        % freedom increases by 1
        function addStretchingMode(this,flag)
            this.stretchingMode = flag;
            if flag == true
                this.ndof = this.ndof + 1;
            end
        end

        %------------------------------------------------------------------
        % Enables the rotation mode. If enables, the number of degrees of
        % freedom increases by 1
        function addRelRotationMode(this,flag)
            this.relRotationMode = flag;
            if flag == true
                this.ndof = this.ndof + 1;
            end
        end

        %------------------------------------------------------------------
        % Initializes the integration points for the element obtaining the
        % coordinates and weights
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
        % Computes the element stiffness matrix, internal force vector and
        % other optional outputs using numerical integration over the
        % element
        % 
        % Outputs:
        %   Ke    - Element stiffness matrix.
        %   Ce    - Element damping matrix.
        %   fi    - Internal force vector.
        %   fe    - External force vector.
        %   dfidu - Derivative of internal force with respect to 
        %           displacement.
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
        % Computes the enrichment interpolation matrix for the given
        % integration point coordinates, reference point and tangential
        % vector
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