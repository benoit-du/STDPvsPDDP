function STDPvsEbPDDP(N,omega,Delta,theta_0_std,kappa_0,kappa_0_std,tau_plus,tau_R,beta,A_plus,n_modes,outputDir,seed,timeInS,dt)

%%% 27-04-23    first revision
%%% Benoit Duchet, University of Oxford

%%% timeInS         = duration of simulation (biological time in seconds)
%%% dt              = simulation time step (in seconds)

%%% N               = number of oscillators
%%% Omega           = average natural frequency of oscillators (rd/s)
%%% Delta           = standard deviation of oscillator frequency distribution (rd/s)

%%% tau_plus        = timescale of long-term potentiation (LTP) in seconds
%%% tau_R           = tau_minus/tau_plus, i.e. ratio of the timescale of long-term depression (LTD) to the timescale of LTP 
%%% A_plus          = magnitude of LTP
%%% beta            = A_minus/A_plus, i.e. ratio of the magnitude of LTD to LTP    
%%% n_modes         = number of Fourier modes in the PDDP approximation

%%% theta_0_std     = standard deviation of the distribution of initial phases 
%%% kappa_0         = mean of the distribution of initial coupling strengths 
%%% kappa_0_std     = standard deviation of the distribution of initial coupling strength

%%% outputDir       = directory where output will be saved
%%% seed            = seed used to initialise the random number generator (set to a fixed value for repeatable results)


%PDDP pars equivalent to STDDP pars
eps         = A_plus/(2*pi/omega);
A_p         = 1;
A_d         = beta;
tau_p       = omega*tau_plus;
tau_d       = tau_R*tau_p;

useLatex(true)

parStr = ['Omega=' num2str(omega) '_Delta=' num2str(Delta) '_k=' num2str(kappa_0) '_k_std=' num2str(kappa_0_std) '_tau_plus=' num2str(tau_plus) '_tau_R=' num2str(tau_R) '_beta=' num2str(beta) '_A_p=' num2str(A_plus)];
dirStr = [outputDir filesep 'STDPvsEbPDDP_' parStr '_' getDateTimeStr];
mkdir(dirStr);

n = round(timeInS/dt);

%%% packing parameters
par.N = N;
par.n_p = 1;
%oscillator distrib
par.omega = omega;
par.Delta = Delta;
%weights
par.k11 = kappa_0;%
par.k11_std = kappa_0_std;
par.kSamplers.k11 = @(sz) normrnd(kappa_0,kappa_0_std,sz);
%initial phase distrib
par.psi_1_0 = 0;
par.theta_0_std = theta_0_std;
%stdp
par.beta = beta;
par.A_plus = A_plus;
par.tau_plus = tau_plus;
par.tau_R = tau_R;
%pddp
par.eps = eps;
par.A_p = A_p;
par.A_d = A_d;
par.tau_p = tau_p;
par.tau_d = tau_d;
%Fourier expansion
par.f = getFourierEqRule(eps,A_p,A_d,tau_p,tau_d,n_modes,false);
% simulation parameters
par.repeatable_seed = seed;
simPar.dt = dt;
simPar.n = n;
simPar.snapFact = simPar.n/200;
simPar.plotOmega = true;
simPar.dirStr = dirStr;

%%% side by side comparison of rules
nPts = 500;
phi = linspace(0,2*pi,nPts);
dK_pddp = eps*(2*pi/omega)*(A_p*exp(-phi/tau_p)-A_d*exp((phi-2*pi)/tau_d));

delta_t_p = linspace(0,2*pi/omega,nPts);
dK_stdp_p = A_plus*exp(-delta_t_p/tau_plus);
dK_stdp_m = -A_plus*beta*exp(-abs(delta_t_p-2*pi/omega)/(tau_plus*tau_R));

figure
hold on
plot(phi,dK_pddp,'displayName','pddp')
plot(phi,dK_stdp_p+dK_stdp_m,'--','displayName','stdp')
plot(phi,feval(par.f,phi)*(2*pi/omega),'displayName',['pddp, ' num2str(n_modes) ' Fourier modes'])
legend('location','best')
xlabel('\phi')
fnames{2} = mySaveasFlex('dimXY',[7,7],'fNameNoNowStr',[dirStr filesep 'ruleComp']);
close

%%% simulation of full network with STDP and with PDDP approximation
[~,r_pop_pddp,k_avg_pddp,omegaVect_pddp,snap_pddp,] = kuramoto_onePop_FourierEbPDDP_fwdSim(par,simPar);
[~,r_pop_stdp,k_avg_stdp,omegaVect_stdp,snap_stdp] = kuramoto_onePop_uSTDP_fwdSim(par,simPar);


%%% plotting
ft = 14;
t = (0:simPar.n-1)*simPar.dt;
l_osc = 0.5;

figure
% phases
subplot(2,1,1)
hold on
plot(t,angle(squeeze(r_pop_stdp)),'linewidth',2*l_osc,'displayName','STDP')
plot(t,angle(squeeze(r_pop_pddp)),'linewidth',2*l_osc,'displayName',['ebPDDP, ' num2str(n_modes) ' Fourier modes'])
ylabel('$\psi$','interpreter','latex')
mainTitleStr = {['$\Omega = $ ' num2str(par.omega,2) ', $\Delta = $ ' num2str(par.Delta,2) ...
    ', $\kappa = $ ' num2str(par.k11) ',  $\sigma_\kappa = $ ' num2str(kappa_0_std)],...
    ['$A_+ = $ ' num2str(A_plus) ', $\beta = $ ' num2str(beta) ', $\tau_+ = $ ' num2str(tau_plus) ...
    ', $\tau_- = $ ' num2str(tau_R*tau_plus)], ['ebPDDP with ' num2str(n_modes) ' Fourier modes']};
title(mainTitleStr,'interpreter','latex')
set(gca,'fontsize',ft)
% amplitudes
rho_stdp = abs(squeeze(r_pop_stdp));
rho_pddp = abs(squeeze(r_pop_pddp));
subplot(2,1,2)
hold on
plot(t,rho_stdp,'linewidth',2*l_osc,'displayName','STDP')
plot(t,rho_pddp,'linewidth',2*l_osc,'displayName',['ebPDDP, ' num2str(n_modes) ' Fourier modes'])
ylim([0 1])
legend('location','best')
xlabel('time (s)','interpreter','latex')
ylabel('$\rho$','interpreter','latex')
set(gca,'fontsize',ft)
fnames{1} = mySaveasFlex('dimXY',[16,14],'fNameNoNowStr',[dirStr filesep 'compPhaseAmp']);
close

% average weights
k_avg_stdp_sq = squeeze(k_avg_stdp);
k_avg_pddp_sq = squeeze(k_avg_pddp);
figure
hold on
plot(t,k_avg_stdp_sq,'linewidth',2*l_osc,'displayName','STDP')
plot(t,k_avg_pddp_sq,'linewidth',2*l_osc,'displayName','ebPDDP')
ylabel('$\hat{\kappa}$','interpreter','latex')
title(['$\Omega = $ ' num2str(par.omega,2) ', $\Delta = $ ' num2str(par.Delta,2) ', seed = ' num2str(seed)],'interpreter','latex')
set(gca,'fontsize',ft)
legend('location','eastoutside')
xlabel('time (s)','interpreter','latex')
pause(1)
fnames{3} = mySaveasFlex('dimXY',[12,10],'fNameNoNowStr',[dirStr filesep 'compWeights']);
close

% coupling distribution over time
ft_k = 17;
[c_max_stdp,~,~,k_min_stdp,k_max_stdp] = plotWeightDistrOverTime(snap_stdp,k_avg_stdp,t,simPar,par,N,dirStr,'STDP',ft_k,false,[],[],[]);
plotWeightDistrOverTime(snap_pddp,k_avg_pddp,t,simPar,par,N,dirStr,['ebPDDP, ' num2str(n_modes) ' Fourier modes'],ft_k,false,[],[],[]);

[~,fnames{4}] = plotWeightDistrOverTime(snap_stdp,k_avg_stdp,t,simPar,par,N,dirStr,'STDP',ft_k,true,c_max_stdp,k_min_stdp,k_max_stdp);
[~,fnames{6}] = plotWeightDistrOverTime(snap_pddp,k_avg_pddp,t,simPar,par,N,dirStr,['ebPDDP, ' num2str(n_modes) ' Fourier modes'],ft_k,true,c_max_stdp,k_min_stdp,k_max_stdp);

% coupling matrices
% sortType: 0 is no sorting, 1 by natural freq, 2 by phase velocity
n_t_plot = 4;%plotting four time points
sortType = 1;
fnames{5} = plotCouplingMat(snap_stdp,N,n_t_plot,sortType,omegaVect_stdp,[],dirStr,'STDP',ft_k,[k_min_stdp k_max_stdp]);
fnames{7} = plotCouplingMat(snap_pddp,N,n_t_plot,sortType,omegaVect_pddp,[],dirStr,['ebPDDP, ' num2str(n_modes) ' Fourier modes'],ft_k,[k_min_stdp k_max_stdp]);

%%% summary figure
toShow = imtile(fnames,'ThumbnailSize',[2500 2500],'BackgroundColor','w');
imshow(toShow);
title(mainTitleStr,'interpreter','latex','fontsize',24)
pause(0.8)
mySaveasFlex('dimXY',[40 40],'fNameNoNowStr',[outputDir filesep parStr]);
close

end