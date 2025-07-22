%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuous time incomplete markets codes (Bewley-Aiyagari-Huggett)
%
% Author: Alex Clymo
% Date: July 2025
% Repository: github->alexclymo->incomplete-markets
%
% Solves the steady state policies and distributions and calibrates
% parameters. Heavily based on the original Achdou et al (2022) codes, with
% some tweaks for quickly building the A matrix for larger problems, and
% options to select model variants. 
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

%choose calibration version (huggett, aiyagari, or pe)
cali_version = 'aiyagari'

%load initial guesses or not
load_init = 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other fixed controls (i.e. not changed often or between calibrations)

%number of nodes
Na = 200; %size of asset grid
        
%calibration tolerance
tolC = 1e-5;

%value function tol (should be stricter than tolC)
tolV = 1e-8;
maxIterV = 1000;

% dampening parameters
Delta = 10; %HJB worker HJB dampening (lower is more dampening)
zetaV = 1;  %extra dampening (lower is more dampening)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed parameters

rho = -log(1-0.05); %discount rate
sigma = 2; %risk aversion coefficient, u(c) = c^(1-sigma)/(1-sigma);
aBarMult = -3/12; %borrowing constraint (aBar < 0 = borrowing) as multiple of average yearly income


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters, targets, and settings for each calibration version

switch cali_version
    case 'huggett'
        zProc = 'urisk'; %choose income process: twostate, urisk, ar1
        wss = 1; %fixed wage per unit of worker productivity
        Ass_supply = 0; % assets in zero net supply
        aMaxMult = 0.5; %maximum asset node as multiple of average yearly income
    case 'aiyagari'
        zProc = 'twostate'; %choose income process: twostate, urisk, ar1
        alpha = 1/3; %capital share
        delta = -log(1-0.1); %capital depreciation
        wss_target = 1;
        aMaxMult = 12; %maximum asset node as multiple of average yearly income
    case 'pe'
        zProc = 'ar1'; %choose income process: twostate, urisk, ar1
        wss = 1;
        rss = -log(1-0.03);
        aMaxMult = 12; %maximum asset node as multiple of average yearly income
    otherwise
        error('invalid cali version')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guesses and setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Productivity grid

switch zProc
    case 'twostate' %arbitary two state process 
        Nz = 2; %two state process
        dz = 0.5;
        zgrid = [1-dz;1+dz];
        alphaZ = 1; %shocks arrive roughly once a year
        g12 = 0.4; g21 = 0.4;
        gammaZ = [1-g12,g12;g21,1-g21];
    case 'urisk' %two state process representing employed/unemployed
        Nz = 2; %two state process
        alphaZ = 12; %shocks arrive roughly once a month
        EU_target = 0.03; %target for monthly EU rate
        U_target = 0.06; %target for unemployment rate
        RR_target = 0.4; %replacement rate while unemployed
        % EU_target = 1 - exp(-EU/12)
        EUy = -12*log(1-EU_target); %EU rate as yearly flow
        % EU*(1-U) = UE*U
        UEy = EUy*(1-U_target)/U_target; %UE rate as yearly flow
        if alphaZ < max([EUy,UEy]), error('alphaZ too low'), end
        % build gammaZ and zgrid. alphaZ*gammaZ(1,2) = UEy and
        % alphaZ*gammaZ(2,1) = EUy
        zgrid = [RR_target;1]; %zgrid(1) is unemployment, 2 is employment
        gammaZ = [1-UEy/alphaZ,UEy/alphaZ;EUy/alphaZ,1-EUy/alphaZ];
    case 'ar1' %ar(1) style process: new shock at rate alphaZ drawn from discretised ar(1)
        Nz = 11; %size of z grid
        alphaZ = 1; %arrival rate of new prod level
        rhoZ = 0.8; %rho of Z once shock arrives
        sigmaZ = 0.3; %sigma of Z once shock arrives
        % z grid, assuming log mean = 0
        [P_Rouw, z_Rouw] = rouwen(rhoZ, 0, sigmaZ, Nz);
        gammaZ = P_Rouw';
        zgrid = exp(z_Rouw);
    otherwise
        error('invalid z process')
end

%find ergodic distribution of z process
[V,D] = eig(gammaZ'); %eigenvalues and vectors V and D
i1 = find(abs(diag(D)-1) < 1e-10,1,'first'); %eigenval = 1
if ~isempty(i1)
    gammaZ_erg = V(:,i1)/sum(V(:,i1)); %ergodic dist, scale sum=1
end

%actual expected value of discretised z:
Ez = zgrid'*gammaZ_erg;
stdz = sqrt( (zgrid'-Ez).^2*gammaZ_erg );
Elogz = log(zgrid)'*gammaZ_erg;
stdlogz = sqrt( (log(zgrid)'-Elogz).^2*gammaZ_erg );

%Lss is equal to Ez
Lss = Ez;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Asset grid: agrid = aBar + xgrid

%aBar = aBarMult * mean labour income
aBar = aBarMult * Lss * Ez;
aMax = aMaxMult * Lss * Ez;
%equally spaced grid for assets
agrid = linspace(aBar, aMax, Na)'; 

%normalised asset grid: xgrid = agrid - aBar (useful for borrowing
%constraint shock transition experiments)
xgrid = agrid - aBar;
xMin = min(xgrid); xMax = max(xgrid);

%forward and backwards differences for asset grid (same for xgrid and agrid)
deltaXF = [xgrid(2:end)-xgrid(1:end-1) ; xgrid(end)-xgrid(end-1)]; %forward differences
deltaXB = [xgrid(2)-xgrid(1) ; xgrid(2:end)-xgrid(1:end-1)]; %backward differences
%trapezoid integration nodes for asset grid (only needed for this grid
%since z is discrete by assumption
trX = 0.5*[deltaXF(1) ; deltaXF(2:end-1)+deltaXB(2:end-1) ; deltaXB(end)]; %trapezoid rule weight vectors for x

%grids replicated onto (x,z,rho)
agrid_s = repmat(agrid,[1,Nz]);
xgrid_s = agrid_s - aBar;
deltaXF_s = repmat(deltaXF,[1,Nz]);
deltaXB_s = repmat(deltaXB,[1,Nz]);
trX_s = repmat(trX,[1,Nz]);
zgrid_s = permute(repmat(zgrid,[1,Na]),[2,1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial parameter guesses

switch cali_version
    case 'huggett'
        % initial guess some interest rate below rho
        rss = rho - 0.02;
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
    case 'pe'
        % nothing to guess
    otherwise
        error('invalid cali version')
end

% initial guesses for worker value functions
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
%the true number after the first iteration, so can put any number here
NelemA = 2;

%set up par structure
par.tolV = tolV; par.maxIterV = maxIterV; par.zetaV = zetaV;
par.NelemA = NelemA; par.zgrid_s = zgrid_s; par.agrid_s = agrid_s;
par.Na = Na; par.Nz = Nz; par.deltaXF_s = deltaXF_s;
par.deltaXB_s = deltaXB_s; par.rowInd_A_Xup = rowInd_A_Xup;
par.colInd_A_Xup = colInd_A_Xup; par.rowInd_A_Xdown = rowInd_A_Xdown;
par.colInd_A_Xdown = colInd_A_Xdown; par.sigma = sigma;
par.alphaZ = alphaZ; par.gammaZ = gammaZ; par.rho = rho;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main algorithm: Loop to jointly calibrate and solve general equilibrium
tic

errC = 1; %errC: joint GE and calibration error
errV = 1; %errV: value function error
iterC = 1;
while errC > tolC || errV > tolV %while loop for GE and calibration

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
    
    [v,c,aDot,A,errV,NelemA] = solveWorkerProblem(v,wss,rss,Delta,par);
    par.NelemA = NelemA; %update number of assigments to A matrix

    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve ergodic distribution
    
    %transpose original A matrix
    AT = A';

    %solve ergodic distribution
    ergDist = solveErgDist(AT,trX_s,Na,Nz);

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
            % note: this can be done more simply by iterating over rss as
            % in the huggett version, imposing wss=1 and reverse
            % engineering the implied Kss and Zss each iteration. i iterate
            % over Kss instead to show a different style of update
            Kss_ = Ass;
            errK = (Kss_-Kss)/(1e-5+abs(Kss))
            errW = (wss-wss_target)/(1e-5+abs(wss_target))
            errC = max(abs([errK;errW]));
            % only update parameters if err greater than tolerance
            if errC > tolC
                % adjust Kss towards new value
                zetaK = 0.1; %0 < zetaK < 1. Lower if convergence problems
                Kss = zetaK * Kss_ + (1-zetaK) * Kss;
                % adjust Zss downwards if wss too high
                zetaW = 0.5; %0 < zetaW. Lower if convergence problems
                Zss = Zss - zetaW*errW;
            end
        case 'huggett'
            errA = Ass-Ass_supply
            errC = abs(errA);
            % only update parameters if err greater than tolerance
            if errC > tolC
                % adjust rss downwards if asset demand too high
                zetaR = 0.1; %0 < zetaR. Lower if convergence problems
                rss = rss - zetaR*errA;
            end
        case 'pe'
            errC = 0;
        otherwise
            error('invalid cali version')
    end

    iterC = iterC+1;

end %end calibration while loop


disp(['time to solve GE and calibrate model: ',num2str(toc),' seconds'])
disp(['mass at top node of asset grid: ',num2str(100*gtil_mx(end)/sum(gtil_mx)),'%'])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make any more variables desired post convergence

switch cali_version
    case 'aiyagari'
        Yss = Zss * Kss^alpha * Lss^(1-alpha);
        Iss = delta*Kss;
        IYss = Iss/Yss;
end


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