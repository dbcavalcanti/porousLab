%% Anl class
%
% This in an abstract class that defines the analysis object
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
classdef Anl < handle

    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        type     = [];      
    end

    %% Constructor method
    methods
        function this = Anl(type)
            if nargin > 0
                this.type   = type;
            end
        end
    end
    methods
        %------------------------------------------------------------------
        % Partition and solve a linear system of equations.
        %  f --> free d.o.f. (natural B.C. - unknown) 
        %  c --> constrained d.o.f. (essential B.C. - known) 
        %
        % [ Kff Kfs ] * [ Uf ] = [ Fext ]
        % [ Ksf Kss ]   [ Us ] = [   R  ]
        %
        function [U,Fext] = solveSystem(~,mdl,K,Fext,U)

            if nargin < 5
                U = zeros(mdl.ndof,1);
            end

            % Partition system of equations
            Kff = K(mdl.doffree, mdl.doffree);
            Ff  = Fext(mdl.doffree);
            
            % Solve the system of equilibrium equations
            Uf = Kff \ Ff;

            % Displacement vector
            U(mdl.doffree)  = Uf;
            
        end

        %------------------------------------------------------------------
        % Get the "force" vector associated with the prescribed degrees of
        % freedom.
        function F = addVectorFromPrescDofs(~, K, F, mdl)

            % Get the free and fixed dofs vector
            freedof  = 1:mdl.ndoffree;
            fixeddof = (1+mdl.ndoffree):mdl.ndof;

            % Get the prescribed dofs
            Us = mdl.U(fixeddof);

            % Get the "sub-stiffness" matrix 
            Kfs = K(freedof,fixeddof);

            % Compute the force vector
            F(freedof) = F(freedof) - Kfs * Us;
            
        end
    end
end