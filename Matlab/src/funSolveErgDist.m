function ergDist = funSolveErgDist(AT,trX_s,Na,Nz)

NA = Na*Nz;

%time step in simul (can be huge)
dt = 1000000;

%initialise gtil: this is the density multiplied by agrid trapeze weights.
%should use this for any integrations and calcuations of aggregates
gtil = zeros(NA,1);
gtil(1) = 1;

%simulate
diff = 1;
tol = 1e-6;
while diff > tol
    %update and diff
    gtil_ = (speye(NA) - dt*AT)\gtil;
    gtil_ = max(gtil_,0); %possible to get slight -ve in matrix inversion
    diff = max(abs(gtil-gtil_));
    gtil = gtil_;
end
if abs(sum(gtil)-1)>1e-5, warning('gtil not equal 1'), end
gtil = gtil/sum(gtil);

%reshape to grid size. 
gtil = reshape(gtil,[Na,Nz]);

%now scale whole dist using trapezoid integration in x dimension to get the
%underlying density (never really need to use this)
g = gtil./trX_s;

%marginal distributions among employed across x, z, and rho
gtil_mx = sum(gtil,2); %multiplied by trap weights
gtil_mz = sum(gtil,1);


%pack into structure
ergDist.gtil = gtil;
ergDist.g = g;
ergDist.gtil_mx = gtil_mx;
ergDist.gtil_mz = gtil_mz;

