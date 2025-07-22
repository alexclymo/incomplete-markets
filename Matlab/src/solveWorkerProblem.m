function [v,c,aDot,A,errV,NelemA] = solveWorkerProblem(vInit,w,r,Delta,par)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract par structure and set up

% initial value
v = vInit;

% extrac par structure
tolV = par.tolV;
maxIterV = par.maxIterV;
zetaV = par.zetaV;
NelemA = par.NelemA;
zgrid_s = par.zgrid_s;
agrid_s = par.agrid_s;
sigma = par.sigma;
Na = par.Na;
Nz = par.Nz;
deltaXF_s = par.deltaXF_s;
deltaXB_s = par.deltaXB_s;
rowInd_A_Xup = par.rowInd_A_Xup;
colInd_A_Xup = par.colInd_A_Xup;
rowInd_A_Xdown = par.rowInd_A_Xdown;
colInd_A_Xdown = par.colInd_A_Xdown;
alphaZ = par.alphaZ;
gammaZ = par.gammaZ;
rho = par.rho;

dvfwd = zeros(Na,Nz); %preallocation, forward differenced value function
dvbck = zeros(Na,Nz); %preallocation, backward differenced value function


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Value function iteration


errV = 1;
iterV = 1;
while errV > tolV && iterV <= maxIterV %inner value function loop

    %restart A construction: A is a Na*Nz square matrix
    %IMPORTANT: nesting order for Na-Nz blocks must follow
    %the matlab flattening control (:). This means A is blocks of a
    %nested in blocks of z

    %initialise sparse matrix building vectors
    storeRowInd = zeros(NelemA,1);
    storeColInd = zeros(NelemA,1);
    storeFlows = zeros(NelemA,1);
    %running position of where to store next nodes
    storePos = 0;
    

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OPTIMAL CONSUMPTION AND ASSET DRIFT

    cStay = w*zgrid_s + r*agrid_s;%consumption when assets not changed
    mUStay = cStay.^(-sigma); %marginal utility when assets not changed
    dvfwd(1:Na-1,:,:) = (v(2:Na,:,:) - v(1:Na-1,:,:))./deltaXF_s(1:Na-1,:,:);
    dvfwd(Na,:,:) = mUStay(Na,:,:);
    dvbck(2:Na,:,:) = (v(2:Na,:,:) - v(1:Na-1,:,:))./deltaXB_s(2:Na,:,:);
    dvbck(1,:,:) = mUStay(1,:,:);
    % optimal c from foc
    cfwd = dvfwd.^(-1/sigma);
    cfwd(Na,:,:) = cStay(Na,:,:); %at top, force xDot=0 if try to move up
    cbck = dvbck.^(-1/sigma);
    cbck(1,:,:) = cStay(1,:,:); %at bottom, force xDot=0 if try to move down
    % correct if imaginary number shows up during convergence
    cfwd(abs(imag(cfwd))>0) = 10;
    cbck(abs(imag(cbck))>0) = 10;
    if max(abs(imag(cfwd)),[],'all')>0, error('imag cfwd'), end
    if max(abs(imag(cbck)),[],'all')>0, error('imag cfwd'), end
    % drift in x at optimal c choice in each direction
    mXfwd = w*zgrid_s + r*agrid_s - cfwd;
    mXbck = w*zgrid_s + r*agrid_s - cbck;
    % check which direction is optimally taken
    iBoth = mXfwd>0 & mXbck < 0; %both valid
    iFwdOnly = mXfwd > 0 & mXbck >=0; %only forward valid
    iBckOnly = mXfwd<=0 & mXbck < 0; %only back valid
    if sigma ~= 1 %forward and back hamiltonian
        hFwd = cfwd.^(1-sigma)/(1-sigma) + dvfwd.*mXfwd;
        hBck = cbck.^(1-sigma)/(1-sigma) + dvbck.*mXbck;
    else
        hFwd = log(cfwd) + dvfwd.*mXfwd;
        hBck = log(cbck) + dvbck.*mXbck;
    end
    indFwd = iFwdOnly + iBoth.*(hFwd>=hBck); %forward is optimal
    indBck = iBckOnly + iBoth.*(hFwd<hBck); %back is optimal
    indStay = 1 - indFwd - indBck; %staying is optimal
    % consumption, utility, and x drift at optimal choice
    c = cfwd.*indFwd + cbck.*indBck + cStay.*indStay;
    xDot = mXfwd.*indFwd + mXbck.*indBck;

    % store asset contribution to diagonal term and off diagonals
    diag_x = -mXfwd./deltaXF_s.*indFwd + mXbck./deltaXB_s.*indBck;
    alphaUp = mXfwd./deltaXF_s.*indFwd;
    alphaDown = -mXbck./deltaXB_s.*indBck;

    % add asset off diagonals to A matrix
    nNodes = length(rowInd_A_Xup); %number of nodes adding to the list
    storeRowInd(storePos + (1:nNodes)) = rowInd_A_Xup;
    storeColInd(storePos + (1:nNodes)) = colInd_A_Xup;
    temp = alphaUp(1:end-1,:);
    storeFlows(storePos + (1:nNodes)) = temp(:);
    storePos = storePos + nNodes;

    nNodes = length(rowInd_A_Xdown); %number of nodes adding to the list
    storeRowInd(storePos + (1:nNodes)) = rowInd_A_Xdown;
    storeColInd(storePos + (1:nNodes)) = colInd_A_Xdown;
    temp = alphaDown(2:end,:);
    storeFlows(storePos + (1:nNodes)) = temp(:);
    storePos = storePos + nNodes;


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add z shock transitions

    %z shock
    for i = 1:Nz
        %for each i, this is the flow to i'
        temp = alphaZ*repmat(gammaZ(i,:),[Na,1]);

        nNodes = length(temp(:)); %number of nodes adding to the list
        storeRowInd(storePos + (1:nNodes)) = repmat( (1:Na)' + Na*(i-1),[Nz,1]);
        storeColInd(storePos + (1:nNodes)) = (1:Na*Nz)';
        storeFlows(storePos + (1:nNodes)) = temp(:);
        storePos = storePos + nNodes;

    end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finish building A matrix

    % U diag on (x)
    diagg = diag_x(:) - alphaZ;

    nNodes = length(diagg); %number of nodes adding to the list
    storeRowInd(storePos + (1:nNodes)) = (1:Na*Nz)';
    storeColInd(storePos + (1:nNodes)) = (1:Na*Nz)';
    storeFlows(storePos + (1:nNodes)) = diagg;
    storePos = storePos + nNodes;


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INVERT A MATRIX AND UPDATE

    %trim storeFlows if NelemA ever happened to shrink (should never
    %happen, but just in case)
    storeRowInd = storeRowInd(1:storePos);
    storeColInd = storeColInd(1:storePos);
    storeFlows = storeFlows(1:storePos);

    % construct A matrix
    A = sparse(storeRowInd, storeColInd, storeFlows, Na*Nz,Na*Nz);

    %update total number of non-empty nodes for storeFlows construction
    %in next iteration
    NelemA = storePos;

    %error check: A should sum to zero along each row
    if any(full(sum(A,2)) > 1e-10), error('A rows don''t sum to 0'), end

    %flow utility of unemployed and employed given optimal consumption
    if sigma ~= 1
        u = c.^(1-sigma)/(1-sigma);
    else
        u = log(c);
    end

    %matrix inversion to update values
    B = (rho + 1/Delta)*speye(Na*Nz) - A;
    uStack = u(:);
    vStack = v(:);

    bStack = uStack + vStack/Delta;
    vNew = B\bStack;

    errV = max(abs((vStack(:) - vNew(:))./vStack(:)));
    iterV = iterV + 1;

    if any(abs(imag(vNew(:)))), error('imaginary numbers in v'), end

    % Update v with Extra dampening if needed
    vNew = reshape(vNew,[Na,Nz]);
    if errV > tolV
        v = zetaV*vNew + (1-zetaV)*v;
    end


end


if iterV == maxIterV + 1
    disp(' ')
    disp(['WARNING: worker vfi hit max iter. errV = ',num2str(errV)])
    disp(' ')
end

aDot = xDot;