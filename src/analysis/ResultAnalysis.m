%% ResultAnalysis Class
%
% A result analysis object is responsible for storing the analysis results.
%
classdef ResultAnalysis < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Result options
        name    = [];  % vector of curves names
        node    = [];  % vector of curves nodes
        dof     = [];  % vector of curves d.o.f's
        coordP  = [];  % Coordinate of the point to plot the time evolution
        coordPf  = []; % Coordinate of the point to plot the time evolution
        coordST = [];  % Coordinate of the point to plot the time evolution
        
        % Result storage
        steps = 0;     % number of performed steps
        niter = 0;     % number of total iterations
        lbd   = [];    % vector of load ratios of all steps
        F     = [];    % vector of nodal forces and reactions
        U     = [];    % matrix of nodal displacement vectors of all steps/modes
        time  = [];    % Vector with the time 
        p     = [];    % Vector with the pressure
        pf    = [];
        ST    = [];    % Vector with the slip tendency
        tn    = [];    % Vector with the shear cohesive stress
        ts    = [];    % Vector with the normal cohesive stress
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function res = ResultAnalysis(dof,coordP, coordPf, coordST)
            res.dof = dof;
            res.coordP  = [];
            res.coordPf = [];
            res.coordST = [];
            if isempty(coordP) == false
                res.coordP = coordP;
            end
            if isempty(coordPf) == false
                res.coordPf = coordPf;
            end
            if isempty(coordST) == false
                res.coordST = coordST;
            end
        end
    end

    %% Public methods
    methods
        function plotCurves(this)
            figure
            hold on
            box on, grid on, axis on
            plot(this.U, this.lbd, 'o-k');
            %plot([0.0,1.0,4.5],[0.0,0.5,0.0],'-b')
            %legend('Numerical','Analytical','Interpreter','latex')
            xlabel('Displacement (mm)','Interpreter','latex')
            ylabel('Load factor','Interpreter','latex')
            xaxisproperties= get(gca, 'XAxis');
            xaxisproperties.TickLabelInterpreter = 'latex'; 
            yaxisproperties= get(gca, 'YAxis');
            yaxisproperties.TickLabelInterpreter = 'latex';   
            set(gca,'FontSize',14);
        end

        function plotTimeEvolution(this)
            figure
            hold on
            box on, grid on, axis on
            % days = 24*60*60;
            days = 1.0;
            plot(this.time/days, this.p, '-k','LineWidth',1.5);
            xlabel('Time (s)')
            ylabel('Pressure (kPa)')
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
        end

        function plotDiscontinuityPressureTimeEvolution(this)
            figure
            hold on
            box on, grid on, axis on
            % days = 24*60*60;
            days = 1.0;
            plot(this.time/days, this.pf, '-k','LineWidth',1.5);
            xlabel('Time (s)')
            ylabel('Discontinuity pressure (kPa)')
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
        end

        function plotSlipTendencyEvolution(this)
            figure
            hold on
            box on, grid on, axis on
            days = 24*60*60;
            for i = 1:size(this.ST,1)
                plot(this.time/days, this.ST(i,:),'LineWidth',1.5);
                legendInfo{i} = ['y = ' num2str(round(this.coordST(i,2),2))];
            end
            xlabel('Time (days)')
            ylabel('Slip tendency')  
            legend(legendInfo);
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
        end

        function plotSlipTendencyVsPressure(this)
            figure
            hold on
            box on, grid on, axis on
            for i = 1:size(this.ST,1)
                plot(this.p/1000, this.ST(i,:),'LineWidth',1.5);
                legendInfo{i} = ['y = ' num2str(round(this.coordST(i,2),2))];
            end
            xlabel('Pressure (MPa)')
            ylabel('Slip tendency')  
            legend(legendInfo);
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
        end

        function plotNormalCohesiveStressEvolution(this)
            figure
            hold on
            box on, grid on, axis on
            days = 24*60*60;
            for i = 1:size(this.tn,1)
                plot(this.time/days, this.tn(i,:),'LineWidth',1.5);
                legendInfo{i} = ['y = ' num2str(round(this.coordST(i,2),2))];
            end
            xlabel('Time (days)')
            ylabel('Normal cohesive stress (kPa)')  
            legend(legendInfo);
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
        end

        function plotShearCohesiveStressEvolution(this)
            figure
            hold on
            box on, grid on, axis on
            days = 24*60*60;
            for i = 1:size(this.ts,1)
                plot(this.time/days, this.ts(i,:),'LineWidth',1.5);
                legendInfo{i} = ['y = ' num2str(round(this.coordST(i,2),2))];
            end
            xlabel('Time (days)')
            ylabel('Shear cohesive stress (kPa)')  
            legend(legendInfo);
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
        end

        function plotCohesiveStressesEvolution(this,c,phi,ft)
            figure
            hold on
            box on, grid on, axis on
            for i = 1:size(this.tn,1)
                plot(this.tn(i,:), this.ts(i,:),'LineWidth',1.5);
                legendInfo{i} = ['y = ' num2str(round(this.coordST(i,2),2))];
            end
            xline(ft,'--r')
            tnres = linspace(min(min(this.tn)), c/tan(phi),100);
            tsres = c + tan(phi) * tnres;
            plot(tnres,  tsres, '--r','LineWidth',1.0);
            plot(tnres, -tsres, '--r','LineWidth',1.0);
            xlabel('Normal cohesive stress (kPa)')
            ylabel('Shear cohesive stress (kPa)')  
            legend(legendInfo);
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
        end

    end
end