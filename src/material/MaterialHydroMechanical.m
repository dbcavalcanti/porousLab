%% Material class
%
% This class defines an abstract stress-strain constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef MaterialHydroMechanical < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        model       = 'elastic';   
        parameters  = [];
        anm         = 'PlaneStress';
        nStVar      = 0;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialHydroMechanical(model, parameters, anm)
            this.model      = model;
            this.parameters = parameters;
            this.anm        = anm;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Get the number of state variables associated with the model 
        % besides the strains and stresses, like plastic deformations or
        % damage
        function nStVar = getNumberStateVar(this)
            nStVar = this.nStVar;
        end

    end

    %% Abstract methods
    methods(Abstract)

        %------------------------------------------------------------------
        % Compute the stress vector
        stress = stressVct(this, strain, pt);

        %------------------------------------------------------------------
        % Compute the constitutive matrix
        De = constitutiveMtrx(this, strain, pt);

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive
        % matrix
        [stress,De] = evalConstitutiveModel(this,strain,pt);
        
    end

    %% Public methods
    methods
        
        %------------------------------------------------------------------
        % Compute the permeability matrix
        function K = permeabilityMtrx(this)

            % Conductivity
            k = this.parameters(3);

            % Water specific weight
            gw = this.waterSpecificWeight();

            % Permeability matrix: K/gw = k/mu
            K = [k 0.0; 0.0 k]/gw;
            % K = [2.0 0.0; 0.0 1.0]*1.0e-8/gw;

        end

        %------------------------------------------------------------------
        % Compute the permeability matrix
        function comp = compressibilityCoeff(this)

            % Get material parameters
            Ks   = this.parameters(6);
            Kf   = this.parameters(7);
            biot = this.parameters(8);
            por  = this.porosity();

            if Ks == 0.0, Ks = 1.0; end
            if Kf == 0.0, Kf = 1.0; end

            % Compressibility coefficient: comp = 1 / M
            comp = (biot - por)/Ks + por/Kf;

        end

        %------------------------------------------------------------------
        % Get the porosity
        function por = porosity(this)

            por = this.parameters(9);

        end

        %------------------------------------------------------------------
        % Get the Biot's coefficient
        function biot = biotCoeff(this)

            biot = this.parameters(8);

        end

        %------------------------------------------------------------------
        % Get the solid specific weight
        function gw = waterSpecificWeight(this)

            gw = max(this.parameters(4),0.00001);

        end

        %------------------------------------------------------------------
        % Get the solid specific mass
        function rhos = solidSpecificMass(this)

            rhos = this.parameters(10);

        end

        %------------------------------------------------------------------
        % Get the solid specific weight
        function gs = solidSpecificWeight(this)

            rhos = this.solidSpecificMass();
            grav = this.gravityAcc();

            gs = rhos * grav / 1000;

        end

        %------------------------------------------------------------------
        % Get the gravity acceleration
        function g = gravityAcc(this)

            g = this.parameters(11);

        end

        %------------------------------------------------------------------
        % Get the bulk specific mass
        function rho = rhoBulk(this)

            % Get material parameters
            por  = this.porosity();
            gw   = this.waterSpecificWeight();
            grav = this.gravityAcc();
            rhos = this.solidSpecificMass();
            
            % Specific mass of the bulk
            % rho =  (por*gw/grav) + rhos/1000.0;
            rho =  rhos/1000.0;

        end


    end
end