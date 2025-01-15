%% Material_Elastic class
%
% This class defines an linear elastic stress-strain constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef RelativePermeabilityUMAT < RelativePermeability  
     %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        Sl_curve = [];          % Must be sorted!!
        kr_curve = [];          % Must be sorted!!
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityUMAT(Sl_curve,kr_curve)
            this = this@RelativePermeability('umat');
            this.Sl_curve = Sl_curve;
            this.kr_curve = kr_curve;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the relative permeability
        function kr = calculate(this, Sl, porousMedia)
            kr = interp1(this.Sl_curve,this.kr_curve,Sl,'linear', 'extrap');
            kr = max(min(kr, 1.0), porousMedia.klrmin);
        end
        
    end
end