%% EnrichedElement_H class
%
% This class defines a single-phase fluid flow finite element 
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef EnrichedElement_H < RegularElement_H   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_H(type, node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric)
            this = this@RegularElement_H(type, node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        function [Ke, Ce, fi, fe, dfidu] = elementData(this)

           % Get the continuum contribution
           [Ke, Ce, fi, fe, dfidu] = elementData@RegularElement_H(this);
           if ~isempty(this.discontinuity)
               [Kd, Sd, Tcd] = getDiscontinuitiesData(this);
               Ke = Ke + Tcd' * Kd * Tcd;
               Ce = Ce + Tcd' * Sd * Tcd;   
           end
        end

        %------------------------------------------------------------------
        function [Kd, Sd, Tcd] = getDiscontinuitiesData(this)

            nEnrDofs          = this.getNumberEnrichedDofs();
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();

            % Initialize the output data 
            Kd = zeros(nEnrDofs,nEnrDofs);
            Sd = zeros(nEnrDofs,nEnrDofs);

            % Condensation matrix
            Tcd = zeros(2*nDiscontinuities,this.nglp);

            % Loop through the discontinuities
            k = 1;
            for i = 1:nDiscontinuities

                % Dofs associated with this discontinuity segment
                dofs = nDofDiscontinuity*(i-1)+1 : nDofDiscontinuity*i;

                % Get the discontinuity data
                [Kdi,Sdi,~,~,~] = this.discontinuity(i).elementData();

                % Assemble the contribution of this discontinuity
                Kd(dofs,dofs) = Kdi;
                Sd(dofs,dofs) = Sdi;

                % Loop through the nodes of the discontinuity to fill the
                % condensation matrix
                for j = 1:2
                    X = this.discontinuity(i).node(j,:);
                    Xn = this.shape.coordCartesianToNatural(this.node,X);
                    Np = this.shape.shapeFncMtrx(Xn);
                    Tcd(k,:) = Np;
                    k = k + 1;
                end
                
            end
        end

        %------------------------------------------------------------------
        function nEnrDof = getNumberEnrichedDofs(this)
            nEnrDof = this.getNumberOfDiscontinuities();
            nEnrDof = nEnrDof * this.getNumberOfDofPerDiscontinuity();
        end

        %------------------------------------------------------------------
        function n = getNumberOfDiscontinuities(this)
            n = size(this.discontinuity,1);
        end

        %------------------------------------------------------------------
        function n = getNumberOfDofPerDiscontinuity(this)
            n = 2;  
        end

        %------------------------------------------------------------------
        function addDiscontinuitySegment(this,dseg)
            this.discontinuity = [this.discontinuity; dseg];
        end
        %------------------------------------------------------------------
    end
end