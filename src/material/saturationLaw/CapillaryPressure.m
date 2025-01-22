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
classdef CapillaryPressure < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id = 'name1';   
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = CapillaryPressure(matModel)
            this.id = matModel;
        end
    end

    %% Abstract methods
    methods(Abstract)

        %------------------------------------------------------------------
        % Liquid saturation degree
        Sw = saturationDegree(this, pc, porousMedia);

        %------------------------------------------------------------------
        % Derivative of the liquid saturation degree wrt to the capillary
        % pressure
        dSwdPc = derivativeSaturationDegree(this, pc, porousMedia);
        
    end

end