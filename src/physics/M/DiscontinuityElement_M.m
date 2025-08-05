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
        % Displacement jump order
        function n = displacementJumpOrder(this)
            n = 0;
            if (this.stretchingMode || this.relRotationMode) 
                n = 1;
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

        %------------------------------------------------------------------
        % Projection matrix
        function P = projectionMatrix(this)
            n = this.normalVector();
            P = [n(1) , 0.0;
                 0.0  , n(2);
                 0.0  , 0.0;
                 n(2) , n(1)];
        end

        %------------------------------------------------------------------
        % Integrate the polynomial stress interpolation of the continuum
        % along the discontinuity
        function S = intPolynomialStressIntp(this, celem)

            % Initialize variables
            jumpOrder     = this.displacementJumpOrder();
            dimPolyStress = celem.shape.dimPolynomialStressInterp();
            S = zeros(dimPolyStress, jumpOrder + 1);

            % Get the discontinuity geometric properties
            Xr = this.referencePoint();
            m  = this.tangentialVector();
            ld = this.ld();

            % Numerical integration
            for i = 1:this.nIntPoints

                % Numerical integration term. The determinant is ld/2.
                c = 0.5 * ld * this.intPoint(i).w * this.t;

                % Cartesian coordinates of the integration point 
                X = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Continuum stress interpolation polynomial
                p = celem.shape.polynomialStress(X);

                if (jumpOrder == 0)
                    S = S + p * c;
                elseif (jumpOrder == 1)
                    % Tangential coordinate
                    s = m' * (X' - Xr');
                    S = S + [p , s*p] * c;
                end
            end
        end

        %------------------------------------------------------------------
        % Get specified field. Fill the coordinate matrix and the field.
        function [X, f] = getField(this,field)
            if strcmp(field,'Sn')
                [X, f] = this.getCohesiveStresses(2); 
            elseif strcmp(field,'St')
                [X, f] = this.getCohesiveStresses(1);
            elseif strcmp(field,'Dn')
                [X, f] = this.getDisplacementJump(2);
            elseif strcmp(field,'Dt')
                [X, f] = this.getDisplacementJump(1);
            end
        end

        %------------------------------------------------------------------
        % Get cohesive stresses component
        function [X, f] = getCohesiveStresses(this,stressId)
            X = zeros(this.nIntPoints,2);
            f = zeros(this.nIntPoints,1);
            for i = 1:this.nIntPoints
                % Cartesian coordinates of the integration point 
                X(i,:) = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);
                % Get cohesive stresses
                f(i) = this.intPoint(i).stress(stressId);
            end
        end

        %------------------------------------------------------------------
        % Get cohesive stresses component
        function [X, f] = getDisplacementJump(this,strainId)
            X = zeros(this.nIntPoints,2);
            f = zeros(this.nIntPoints,1);
            for i = 1:this.nIntPoints
                % Cartesian coordinates of the integration point 
                X(i,:) = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);
                % Get strain component
                f(i) = this.intPoint(i).strain(strainId);
            end
            
        end


    end
end