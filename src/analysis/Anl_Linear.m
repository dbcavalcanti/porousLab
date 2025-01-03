%% Anl_Nonlinear Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anl to deal with linear-elastic analysis.
%
classdef Anl_Linear < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_Linear(result)
            this = this@Anl('Linear',result);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function process(~,mdl)

            % Add contribution of the nodal forces to the external force
            % vector
            F = mdl.addNodalLoad(mdl.F);     
            
            % Add contribution of the load applied at the middle-plane of
            % the discontinuity
            % F = mdl.addNodalLoadEnrichMidPlane(F);

            % Compute the global stiffness matrix
%             U0 = mdl.U;       % To check the computation of the stresses
            [K,f] = mdl.setLinearSystem(mdl.U,F); 

            % Compute pre-conditioner matrix
            % S = Anl.SIPICpreconditioner(K);

            % Solve linear system
            % mdl.U(mdl.totFreeDof) = S'*((S*K*S')\(S*f));
            mdl.U(mdl.totFreeDof) = K\f;

            % Save final result
            for el = 1:mdl.nelem
                gle = mdl.element(el).type.gle;
                mdl.element(el).type.ue = mdl.U(gle);
            end

            mdl.setLinearSystem(mdl.U,F); 
            mdl.updateStateVar();


        end

    end
    
end