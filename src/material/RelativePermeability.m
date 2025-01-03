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
classdef RelativePermeability < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id = 'name1';   
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeability(matModel)
            this.id = matModel;
        end
    end

    %% Abstract methods
    methods(Abstract)

        %------------------------------------------------------------------
        % Liquid relative permeability
        krL = liquidRelativePermeability(this, Sw, porousMedia);

        %------------------------------------------------------------------
        % Gas relative permeability
        krG = gasRelativePermeability(this, Sw, porousMedia);
        
    end

end