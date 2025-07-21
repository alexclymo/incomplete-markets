%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuous time incomplete markets codes (Bewley-Aiyagari-Huggett)
%
% Author: Alex Clymo
% Date: July 2025
% Repository: github->alexclymo->incomplete-markets
%
% Solves the steady state policies and distributions and calibrates
% parameters. Heavily based on the original Achdou et al (2022) codes, 
% adapted for my own style and workflow. Much work was developed over many
% projects in collaboration with Piotr Denderski. 
%
% Use: Use cali_version to select whether solving partial equilibrium,
% Aiyagari, or Huggett style model. Code will solve steady state and
% calibrate certain parameters, and produce plots for the policy functions
% and distribution. 
% 
% One unit of time = 1 year
%
% Partial equilibrium details: (cali_version = 'pe')
%   Fix rss and wss to arbitrary values
% Huggett details: (cali_version = 'huggett')
%   Fix wss, solve rss which clears the bond market
% Aiyagari details: (cali_version = 'aiyagari')
%   Aggregate production function: Y = Z K^alpha L^(1-alpha)
%   Solve Kss which clears the market, calibrate Z to hit wss = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Housekeeping

clear
close all
clc

% add functions / scripts subfolder to path
addpath('./src')

% This script changes all text interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

% Change default plot linewidth
set(groot, 'DefaultLineLineWidth', 2);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls

%choose calibration version (aiyagari, huggett, or pe)
cali_version = 'aiyagari'

%load initial guesses or not
load_init = 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other fixed controls (i.e. not changed often or between calibrations)

%number of nodes
Na = 200; %size of asset grid
Nz = 3; %size of z grid
        
%calibration tolerance
tolC = 1e-5;

%value function tol (should be stricter than tolC)
tolV = 1e-6;
maxIterV = 1000;

% dampening parameters
Delta = 10; %HJB worker HJB dampening (lower is more dampening)
zetaV = 1;  %extra dampening (lower is more dampening)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed parameters

rho = -log(1-0.05); %discount rate
sigma = 2; %risk aversion coefficient, u(c) = c^(1-sigma)/(1-sigma);
aBarMult = -3/12; %borrowing constraint (aBar < 0 = borrowing) as multiple of av yearly wage
alphaZ = 1; %arrival rate of prod shock
rhoZ = 0.7;
sigmaZ = 0.3;

% x grid limits: agrid = aBar + xgrid
xMin = 0; %should be zero by definition
xMax = 5; %set high enough not to bind in ergodic dist


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters, targets, and settings for each calibration version

switch cali_version
    case 'aiyagari'
        alpha = 1/3; %capital share
        delta = 0.1; %capital depreciation
        wss_target = 1;
    case 'huggett'
        wss = 1; %fixed wage
        Ass_supply = 0; % assets in zero net supply
    case 'pe'
        wss = 1;
        rss = -log(1-0.02);
    otherwise
        error('invalid cali version')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guesses and setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Productivity grid

% z grid, initially assuming log mean = 0
[P_Rouw, z_Rouw] = rouwen(rhoZ, 0, sigmaZ, Nz);
gammaZ = P_Rouw';
zgrid = exp(z_Rouw);
%find ergodic distribution of z process
[V,D] = eig(gammaZ'); %eigenvalues and vectors V and D
i1 = find(abs(diag(D)-1) < 1e-10,1,'first'); %eigenval = 1
if ~isempty(i1)
    gammaZ_erg = V(:,i1)/sum(V(:,i1)); %ergodic dist, scale sum=1
end
%actual expected value of discretised z:
Ez = zgrid'*gammaZ_erg;

%Lss is equal to Ez
Lss = Ez;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Asset grid: agrid = aBar + xgrid

%aBar = aBarMult * mean labour income
aBar = aBarMult * Lss * Ez;

%normalise asset grid: xgrid = agrid - aBar
xgrid = linspace(xMin, xMax, Na)'; %equally spaced grid for excess assets

%implied grid for actual assets
agrid = aBar + xgrid;

%forward and backwards differences
deltaXF = [xgrid(2:end)-xgrid(1:end-1) ; xgrid(end)-xgrid(end-1)]; %forward differences
deltaXB = [xgrid(2)-xgrid(1) ; xgrid(2:end)-xgrid(1:end-1)]; %backward differences
%trapezoid integration nodes for asset grid (only needed for this grid
%since z is discrete by assumption
trX = 0.5*[deltaXF(1) ; deltaXF(2:end-1)+deltaXB(2:end-1) ; deltaXB(end)]; %trapezoid rule weight vectors for x

%grids replicated onto (x,z,rho)
xgrid_s = repmat(xgrid,[1,Nz]);
agrid_s = aBar + xgrid_s;
deltaXF_s = repmat(deltaXF,[1,Nz]);
deltaXB_s = repmat(deltaXB,[1,Nz]);
trX_s = repmat(trX,[1,Nz]);
zgrid_s = permute(repmat(zgrid,[1,Na]),[2,1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial parameter guesses

switch cali_version
    case 'aiyagari'
        % We will normalise wss = wss_target by scaling Z
        wss = wss_target;
        % initial guess some interest rate below rho
        rss = rho - 0.02;
        % implied Kss and Zss for this rss guess
        % wage: w = (1-alpha) Z K^alpha L^(-alpha)
        % interest rate: r = alpha Z K^(alpha-1) L^(1-alpha)
        % w/r = (1-alpha)/alpha * K / L
        Kss = wss/rss * alpha/(1-alpha) * Lss;
        Zss = wss / ( (1-alpha) * Kss^alpha * Lss^(-alpha) );
    case 'huggett'
        % initial guess some interest rate below rho
        rss = rho - 0.02;
    case 'pe'
        
    otherwise
        error('invalid cali version')
end

% initial guesses for worker value funs
if sigma ~= 1
    v = ( (wss*zgrid_s + rss*agrid_s).^(1-sigma) )/(1-sigma)/rho;
else
    v = log( wss*zgrid_s + rss*agrid_s )/rho;
end

%total number of nodes
NA = Na*Nz;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load saved initial guesses

switch load_init
    case 1
        load(['./mat/init_cali_',cali_version])
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-build some A matrix indices, and any final pre allocation

% pre build some A matrix indices
[rowInd_A_Xup,colInd_A_Xup,rowInd_A_Xdown,colInd_A_Xdown] = makeAinds(Na,Nz);

%Number of non-empty nodes in the sparse matrix -- will be overwritten with
%the true number after the first iteration
NelemA = 2;

%set up par structure
par.tolV = tolV;
par.maxIterV = maxIterV;
par.zetaV = zetaV;
par.NelemA = NelemA;
par.zgrid_s = zgrid_s;
par.agrid_s = agrid_s;
par.Na = Na;
par.Nz = Nz;
par.deltaXF_s = deltaXF_s;
par.deltaXB_s = deltaXB_s;
par.rowInd_A_Xup = rowInd_A_Xup;
par.colInd_A_Xup = colInd_A_Xup;
par.rowInd_A_Xdown = rowInd_A_Xdown;
par.colInd_A_Xdown = colInd_A_Xdown;
par.sigma = sigma;
par.alphaZ = alphaZ;
par.gammaZ = gammaZ;
par.rho = rho;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main algorithm: Loop to jointly calibrate and solve general equilibrium
tic

errC = 1; errV = 1;
iterC = 1;
while errC > tolC %outer wage loop

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For Aiyagari model, update r and w to be consistent with Kss guess

    switch cali_version
        case 'aiyagari'
            % implied equilibrium interest rate
            rss = alpha * Zss * Kss^(alpha-1) * Lss^(1-alpha) - delta;
            % implied equilibrium wage
            wss = (1-alpha) * Zss * Kss^alpha * Lss^(-alpha);
    end
    

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve worker problem
    
    [v,c,aDot,A,NelemA] = solveWorkerProblem(v,wss,rss,Delta,par);
    par.NelemA = NelemA;

    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve ergodic distribution
    
    %transpose original A matrix
    AT = A';

    %solve ergodic distribution
    ergDist = funSolveErgDist(AT,trX_s,Na,Nz);

    %unpack ergDist structure
    gtil = ergDist.gtil; %distribution pre-scaled with trapezoid weights
    gtil_mx = ergDist.gtil_mx; %marginal distributions
    gtil_mz = ergDist.gtil_mz;


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute aggregates needed for calibration or GE
    
    % steady state assets = sum (average) across agents
    Ass = agrid'*gtil_mx;

    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Errors and updates
    
    switch cali_version
        case 'aiyagari'
            Kss_ = Ass;
            errK = (Kss_-Kss)/(1e-5+abs(Kss))
            errW = (wss-wss_target)/(1e-5+abs(wss_target))
            errC = max(abs([errK;errW]));
            % only update parameters if err greater than tolerance
            if errC > tolC
                % adjust Kss towards new value
                zetaK = 0.5; %above 0.5 starts to cause problems
                Kss = zetaK * Kss_ + (1-zetaK) * Kss;
                % adjust Zss to move wss towards targeted value
                zetaW = 0.5;
                Zss = Zss - zetaW*errW;
            end
        case 'huggett'
            errA = Ass-Ass_supply
            errC = abs(errA);
            % adjust rss towards new value
            zetaR = 0.1; %above 0.1 starts to cause problems
            rss = rss - zetaR*errA;
        case 'pe'
            errC = 0;
        otherwise
            error('invalid cali version')
    end

    iterC = iterC+1;

end %end calibration while loop


disp(['time to calibrate model: ',num2str(toc),' seconds'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot steady state values and policies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot spy of A matrix

spy(A)
% custom lines to show each block
spyXlist = [0,Na,Na + Na*(1:Nz)];
xticks(spyXlist+0.5)
xticklabels(spyXlist)
yticks(spyXlist+0.5)
yticklabels(spyXlist)
grid on
title('$A$ matrix')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot distributions

figure('position',[100,300,1000,300])
MM = 1; NN = 3;
xlims = [min(agrid),max(agrid)];

subplot(MM,NN,1)
plot(agrid,gtil_mx)
grid on
box
xlabel('assets $a$')
xlim(xlims)
title('Marginal asset distribution')

subplot(MM,NN,2)
plot(zgrid,gtil_mz,'-x')
grid on
box
xlabel('Productivity $z$')
%xlim(xlims)
title('Marginal productivity distribution')

subplot(MM,NN,3)
plot(agrid,gtil)
grid on
box
xlabel('assets $a$')
xlim(xlims)
title('Joint distribution $\mu(a,z)$')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot value and policy functions

figure('position',[100,300,1000,300])
MM = 1; NN = 3;
xlims = [min(agrid),max(agrid)];

subplot(MM,NN,1)
plot(agrid,v)
grid on
box
xlabel('assets $a$')
xlim(xlims)
title('value $v(a,z)$')

subplot(MM,NN,2)
plot(agrid,c)
grid on
box
xlabel('assets $a$')
xlim(xlims)
title('consumption $c(a,z)$')

subplot(MM,NN,3)
plot(agrid,aDot)
grid on
box
xlabel('assets $a$')
xlim(xlims)
title('asset drift $\dot a(a,z)$')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save workspace and initial values

% save initial values for calibration loop
switch cali_version
    case 'aiyagari'
        save(['./mat/init_cali_',cali_version], 'v','Kss','Zss')
    case 'huggett'
        save(['./mat/init_cali_',cali_version], 'v','rss')
    case 'pe'
        save(['./mat/init_cali_',cali_version], 'v')
end

% save workspace
save(['./mat/ws_cali_',cali_version])