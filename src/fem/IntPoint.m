%% IntPoint class
%
% This class defines an integration point object
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: January 2023
%%%
%
%% Class definition
classdef IntPoint < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        X                = [];  % Coordinates of the integration point in the natural coordinate system
        w                = 0.0; % Weight associated to the integration point
        strain           = [];  % Current strain vector
        stress           = [];  % Current stress vector
        plasticstrain    = [];  % Current plastic strain vector
        statevar         = [];  % Current state variables vector
        strainOld        = [];  % Previous strain vector
        stressOld        = [];  % Previous stress vector  
        plasticstrainOld = [];  % Current plastic strain vector
        statevarOld      = [];  % Previous state variables vector
        constitutiveMdl  = [];  % Constitutive model object
        anm              = '';  % Analysis model tag
        nVar             = 4;   % Dimension of the stress and strain vectors
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = IntPoint(X,w,constitutiveMdl)
            if nargin > 0
                this.X                = X;
                this.w                = w;
                this.constitutiveMdl  = constitutiveMdl;
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        %  Initialize analysis model (mechanical part)
        function initializeMechanicalAnalysisModel(this,anm)
            this.anm = anm;
            nStvar   = this.constitutiveMdl.getNumberStateVar();
            if strcmp(anm,'PlaneStress')
                this.nVar = 4;
            elseif strcmp(anm,'PlaneStrain')
                this.nVar = 4;
            elseif strcmp(anm,'AxisSymmetrical')
                this.nVar = 4;
            elseif strcmp(anm,'Interface')
                this.nVar = 2;
            end
            this.strain      = zeros(this.nVar,1);
            this.stress      = zeros(this.nVar,1);
            this.statevar    = zeros(nStvar,1);
            this.strainOld   = zeros(this.nVar,1);
            this.stressOld   = zeros(this.nVar,1);
            this.statevarOld = zeros(nStvar,1);
        end

        %------------------------------------------------------------------
        %  Update the current strain vector
        function updateStrainVct(this)
            this.strainOld = this.strain;
        end

        %------------------------------------------------------------------
        %  Update the current state variable vector
        function updateStateVar(this)
            this.statevarOld = this.statevar;
        end

        %------------------------------------------------------------------
        %  Update the current stress vector
        function updateStressVct(this)
            this.stressOld = this.stress;
        end

        %------------------------------------------------------------------
        %  Get the current constitutive matrix
        function D = getConstitutiveMtrx(this,dStrain)
            D = this.constitutiveMdl.constitutiveMtrx(dStrain,this);
        end

        %------------------------------------------------------------------
        %  Get the current constitutive matrix
        function [stress,D] = mechanicalLaw(this)
            [stress,D] = this.constitutiveMdl.mechanicalLaw(this);
            this.stress = stress;
        end

    end
end