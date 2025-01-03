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
        statevar         = [];  % Current state variables vector
        statevarOld      = [];  % Previous state variables vector
        constitutiveMdl  = [];  % Constitutive model object
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = IntPoint(X,w,constitutiveMdl)
            if nargin > 0
                this.X                = X;
                this.w                = w;
                this.statevar         = [];
                this.constitutiveMdl  = constitutiveMdl;
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        %  Update the current state variable vector
        function updateStateVar(this)
            this.statevarOld = this.statevar;
        end

        %------------------------------------------------------------------
        %  Get the current constitutive matrix
        function D = getConstitutiveMtrx(this,dStrain)
            D = this.constitutiveMdl.constitutiveMtrx(dStrain,this);
        end

        %------------------------------------------------------------------
        %  Get the current constitutive matrix
        function [stress,D] = constitutiveModel(this,strain)
            [stress,D] = this.constitutiveMdl.evalConstitutiveModel(strain,this);
            this.strainVct(strain);
            this.stressVct(stress);
        end

    end
end