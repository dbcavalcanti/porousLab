%% RegularElement_M class
%
% This class defines a mechanical finite element 
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef RegularElement_M < RegularElement    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        glu        = [];            % Displacement dofs
        nglu       = 0;             % Number of regular u-dof
        anm        = 'PlaneStrain'; % Analysis model
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement_M(type, node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress)
            this = this@RegularElement(type, node, elem, t, ...
                mat, intOrder, massLumping, lumpStrategy, ...
                isAxisSymmetric);
            this.glu      = glu;
            this.gle      = glu;
            this.nglu     = length(this.glu);
            this.ngle     = length(this.gle);
            if isPlaneStress
                this.anm = 'PlaneStress';
            end
            if isAxisSymmetric
                this.anm = 'AxisSymmetrical';
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
                constModel = Material_M(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
                intPts(i).initializeMechanicalAnalysisModel(this.anm);
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

            % Initialize the matrices
            Ke = zeros(this.nglu, this.nglu);
            Ce = zeros(this.nglu, this.nglu);

            % Initialize external force vector
            fe = zeros(this.nglu, 1);

            % Initialize the internal force vector
            fi = zeros(this.nglu, 1);
            
            % Vector of the nodal dofs
            u  = this.getNodalDisplacement();

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints
               
                % Compute the B matrix at the int. point and the detJ
                [dNdx, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.shape.BMatrix(dNdx);

                % Compute the strain vector
                this.intPoint(i).strain = Bu * u;

                % Compute the stress vector and the constitutive matrix
                [stress,Duu] = this.intPoint(i).mechanicalLaw();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
                
                % Compute the stiffness sub-matrix
                Ke = Ke + Bu' * Duu * Bu * c;

                % Internal force vector
                fi = fi + Bu' * stress * c;
                
                % Compute the gravity forces
                if (this.mat.porousMedia.gravityOn)
                    fe = this.addGravityForces(fe,this.intPoint(i).X,c);
                end
            end
            
        end

        %------------------------------------------------------------------
        % Add contribution of the gravity forces to the external force vct
        function fe = addGravityForces(this, fe, Xn, c)

            % Get gravity vector
            grav = this.mat.porousMedia.g * this.mat.porousMedia.b;

            % Shape function matrix
            N  = this.shape.shapeFncMtrx(Xn);
            Nu = this.shape.NuMtrx(N);

            % Get the porous matrix density
            rhos = this.mat.porousMedia.getDensity();

            % Compute the contribution of the gravitational forces
            fe = fe + Nu' * rhos * grav * c;
            
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the displacement
        function u = getNodalDisplacement(this)
            u = this.ue(1:this.nglu);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the liquid pressure
        function pl = getNodalPressure(this)
            a = this.nglu + 1;
            b = this.nglu + this.nglp;
            pl = this.ue(a:b);
        end

        %------------------------------------------------------------------
        % Function to compute the displacement field in the element.
        function u = displacementField(this,X)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   u   : displacement vector evaluated in "X"
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);
            Nu = this.shape.NuMtrx(Nm);

            % Displacement dof vector
            uv  = this.getNodalDisplacement();
            
            % Regular displacement field
            u = Nu*uv;
        
        end

        %------------------------------------------------------------------
        % Function to compute the stress field inside a given element
        function stress = stressField(this,X,ue)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);

            %Extrapolation matrix
            Q = ones(3,this.nIntPoints);
            for i = 1:this.nIntPoints
                Q(2,i) = this.intPoint(i).X(1);
                Q(3,i) = this.intPoint(i).X(2);
            end
            P = zeros(3,3);
            for i = 1:this.nIntPoints
                P(1,1) = P(1,1) + 1.0;
                P(1,2) = P(1,2) + this.intPoint(i).X(1);
                P(1,3) = P(1,3) + this.intPoint(i).X(2);
                P(2,1) = P(2,1) + this.intPoint(i).X(1);
                P(2,2) = P(2,2) + this.intPoint(i).X(1) * this.intPoint(i).X(1);
                P(2,3) = P(2,3) + this.intPoint(i).X(1) * this.intPoint(i).X(2); 
                P(3,1) = P(3,1) + this.intPoint(i).X(2);
                P(3,2) = P(3,2) + this.intPoint(i).X(2) * this.intPoint(i).X(1);
                P(3,3) = P(3,3) + this.intPoint(i).X(2) * this.intPoint(i).X(2);
            end
            S = P\Q;

            % Matrix with the stress at the integration points
            % Each column corresponds to a stress component:
            % sxx, syy, and tauxy
            stressIP = zeros(this.nIntPoints,3);
            for i = 1:this.nIntPoints
                stressIP(i,1) = this.intPoint(i).stress(1);
                stressIP(i,2) = this.intPoint(i).stress(2);
                stressIP(i,3) = this.intPoint(i).stress(3);
            end

            % Coefficients for the polynomial approximation 
            c = S * stressIP;

            % Interpolated stress field at the given node
            stress = c' * [1.0 ; Xn(1); Xn(2)];

        end

        %------------------------------------------------------------------
        function sn = stressCylindrical(~,stress,X)

            % Get the stress tensor components
            sx = stress(1);
            sy = stress(2);
            tauxy = stress(3);

            % Compute the angle theta
            theta = atan2(X(2), X(1)); % Angle in radians
            
            % Transform the stresses
            sr = sx * cos(theta)^2 + sy * sin(theta)^2 + 2 * tauxy * cos(theta) * sin(theta);
            stheta = sx * sin(theta)^2 + sy * cos(theta)^2 - 2 * tauxy * cos(theta) * sin(theta);
            taurtheta = (sy - sx) * cos(theta) * sin(theta) + tauxy * (cos(theta)^2 - sin(theta)^2);

            sn = [sr, stheta, taurtheta];
        end

        %------------------------------------------------------------------
        % Function to compute the principal stresses
        function [s1,s2] = principalStress(~,stress)

            % Get the stress tensor components
            sxx = stress(1);
            syy = stress(2);
            sxy = stress(3);

            % Mohr's circle center
            c = (sxx + syy) / 2.0;

            % Mohr's circle radius
            r = sqrt(((sxx - syy)/2.0)^2 + sxy^2);

            % Principal stresses
            s1 = c + r;
            s2 = c - r;

        end
    end
end