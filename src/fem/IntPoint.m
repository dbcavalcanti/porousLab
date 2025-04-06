%% IntPoint Class
% This class defines an integration point object used in finite element
% analysis. It stores the properties and methods required to manage the
% state of an integration point, including its coordinates, weights,
% strain, stress, plastic strain, state variables, and constitutive model.
%
%% Methods
% * *IntPoint*: Constructor to initialize the integration point with 
%               coordinates, weight, and constitutive model.
% * *initializeMechanicalAnalysisModel*: Initializes the mechanical 
%                                        analysis model and allocates 
%                                        memory for strain, stress, and 
%                                        state variables based on the 
%                                        analysis type.
% * *updateStrainVct*: Updates the current strain vector and plastic 
%                      strain vector to their previous states.
% * *updateStateVar*: Updates the current state variable vector to 
%                     its previous state.
% * *updateStressVct*: Updates the current stress vector to its 
%                      previous state.
% * *getConstitutiveMtrx*: Retrieves the current constitutive matrix based 
%                          on the given strain increment.
% * *mechanicalLaw*: Computes the stress and constitutive matrix using 
%                    the constitutive model and updates the current 
%                    stress vector.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef IntPoint < handle    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        X                = [];   % Coordinates of the integration point in the natural coordinate system
        w                = 0.0;  % Weight associated to the integration point
        strain           = [];   % Current strain vector
        stress           = [];   % Current stress vector
        plasticstrain    = [];   % Current plastic strain vector
        statevar         = [];   % Current state variables vector
        strainOld        = [];   % Previous strain vector
        stressOld        = [];   % Previous stress vector  
        plasticstrainOld = [];   % Current plastic strain vector
        statevarOld      = [];   % Previous state variables vector
        constitutiveMdl  = [];   % Constitutive model object
        anm              = '';   % Analysis model tag
        nVar             = 4;    % Dimension of the stress and strain vectors
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = IntPoint(X,w,constitutiveMdl)
            if nargin > 0
                this.X = X;
                this.w = w;
                this.constitutiveMdl = constitutiveMdl;
            end
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Initialize analysis model (mechanical part)
        function initializeMechanicalAnalysisModel(this,anm)
            this.anm = anm;

            if strcmp(anm,'PlaneStress')
                this.nVar = 4;
            elseif strcmp(anm,'PlaneStrain')
                this.nVar = 4;
            elseif strcmp(anm,'AxisSymmetrical')
                this.nVar = 4;
            elseif strcmp(anm,'Interface')
                this.nVar = 2;
            end

            nStvar = this.constitutiveMdl.getNumberStateVar();
            this.strain      = zeros(this.nVar, 1);
            this.stress      = zeros(this.nVar, 1);
            this.statevar    = zeros(nStvar,    1);
            this.strainOld   = zeros(this.nVar, 1);
            this.stressOld   = zeros(this.nVar, 1);
            this.statevarOld = zeros(nStvar,    1);

            if this.constitutiveMdl.hasPlasticStrain()
                this.plasticstrain    = zeros(this.nVar,1);
                this.plasticstrainOld = zeros(this.nVar,1);
            end
        end

        %------------------------------------------------------------------
        % Update current strain vector
        function updateStrainVct(this)
            this.strainOld        = this.strain;
            this.plasticstrainOld = this.plasticstrain;
        end

        %------------------------------------------------------------------
        % Update current state variable vector
        function updateStateVar(this)
            this.statevarOld = this.statevar;
        end

        %------------------------------------------------------------------
        % Update current stress vector
        function updateStressVct(this)
            this.stressOld = this.stress;
        end

        %------------------------------------------------------------------
        % Get current constitutive matrix
        function D = getConstitutiveMtrx(this,dStrain)
            D = this.constitutiveMdl.constitutiveMtrx(dStrain,this);
        end

        %------------------------------------------------------------------
        % Get current constitutive matrix
        function [stress,D] = mechanicalLaw(this)
            [stress,D] = this.constitutiveMdl.mechanicalLaw(this);
            this.stress = stress;
        end
    end
end
