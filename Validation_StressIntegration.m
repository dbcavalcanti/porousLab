%% ===================== Elastic plate problem ============================
%
% Elastic traction of a elastic plate validation problem
%
% Author: Danilo Cavalcanti
%
%% ========================================================================
%
% Initialize workspace
clear
initWorkspace; 

%% ============================= MATERIAL =================================

% Create the porous media
rock = PorousMedia('rock');
rock.mechanical = 'vonMises'; % Elastoplastic with von Mises criteria 
rock.Young = 2.0e8;           % Young modulus (Pa)
rock.nu    = 0.3;             % Poisson ratio
rock.sy0   = 1.0e5;           % Initial yield stress (Pa)
rock.Kp    = 0.0;             % Plastic modulus (Pa)

% Material parameters vector
mat  = struct('porousMedia',rock);

%% ============== INITIALIZE THE INTEGRATION POINT ========================

ip = IntPoint([0.0,0.0],1.0,  Material_M(mat));
ip.initializeMechanicalAnalysisModel('PlaneStrain');

%% ========================= STRAIN INCREMENT =============================

% Direction of the strain increment
dstrain0 = [1.0;  % exx
            0.0;  % eyy
            0.0;  % ezz
            0.0]; % gxy

mag = 1.0e-5;
ninc = 150;

idplot       = 1;
stressplot   = zeros(ninc+1,1);
strainplot   = zeros(ninc+1,1);
pstrainplot  = zeros(ninc+1,1);

for i = 1:ninc
    ip.strain = ip.strainOld + mag * dstrain0;
    ip.stress = ip.mechanicalLaw();
    ip.stressOld = ip.stress;
    ip.strainOld = ip.strain;

    % Store
    stressplot(i+1)  = ip.stress(idplot);
    strainplot(i+1)  = ip.strain(idplot);
    pstrainplot(i+1) = ip.plasticstrain(idplot);
end

%  --- Tensão vs. deformação total ----------------------------------------
figure
grid on, box on, hold on
plot(strainplot,stressplot/rock.sy0,'b-','LineWidth',1.5);
xlabel('$\gamma_{xy}$'); ylabel('$\tau_{xy}$/$\tau_Y$');
set(gca,'fontsize',18,'TickLabelInterpreter','latex');

% --- Tensão vs. deformação plástica --------------------------------------
figure
grid on, box on, hold on
plot(pstrainplot,stressplot/rock.sy0,'r-','LineWidth',1.5);
xlabel('$\gamma_{xy}^p$'); ylabel('$\tau_{xy}$/$\tau_Y$'); 
set(gca,'fontsize',18,'TickLabelInterpreter','latex');

