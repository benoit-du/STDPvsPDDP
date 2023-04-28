function [theta,r_pop,k_avg,omegaVect,snap,k0] = kuramoto_onePop_uSTDP_fwdSim(par,simPar)

%%% 27-04-23    first revision
%%% Benoit Duchet, University of Oxford

%%% simulation parameters
dt = simPar.dt;
n = simPar.n;
minJump = 0.9*2*pi;
dirStr = simPar.dirStr;
snapFact = simPar.snapFact;

%%% model parameters
n_p = par.n_p;
N = par.N;
omega = par.omega;
Delta = par.Delta;
kSamplers = par.kSamplers;
A_plus = par.A_plus;
beta = par.beta;
tau_plus = par.tau_plus;
tau_R = par.tau_R;

%%% initialising theta matrix
theta = NaN(n_p,N,n);% n_pop x n_oscillatorPerPop x n_steps

%%% Von Mises distribution
vm_std = par.theta_0_std;
vm_kappa = 1/vm_std^2;
rng(par.repeatable_seed);
theta(1,:,1) = wrapTo2Pi(vmrand(par.psi_1_0,vm_kappa,[1 N]));

k = NaN(n_p,N,n_p,N);
sz = [1,N,1,N];

for i_p1 = 1:n_p
    for i_p2 = 1:n_p
        k(i_p1,:,i_p2,:) = (i_p1 == 1 & i_p2 == 1)*kSamplers.k11(sz);
    end
end

k0 = k;

k_avg = NaN(n_p,n_p,n);% pre x post x time
k_avg(:,:,1) = mean(k,[2,4]);

%%% last spiking index matrix
last_spike_idx = -1E6*ones(n_p,N);

%%% drawing oscillator frequencies
omegaVect = normrnd(omega,Delta,[n_p N]);

if simPar.plotOmega
    nDelta = 5;
    nBins = 15;
    figure
    for i_p = 1:n_p
        subplot(1,n_p,i_p)
        edges = linspace(max(0,omega-nDelta*Delta),omega+nDelta*Delta,nBins);
        histogram(omegaVect(i_p,:),edges)
        xlabel('$\omega$','interpreter','latex')
        title(['pop ' num2str(i_p) ', N = ' num2str(N)],'interpreter','latex')
        set(gca,'fontsize',14)
    end
    mySaveasFlex('dimXY',[12,6],'fNameNoNowStr',[dirStr filesep 'omega']);
    close
end

%%% initialising snapshot matrices
n_snap = round(n/snapFact);
snap.theta = NaN(n_p,N,n_snap);
snap.k = NaN(n_p,N,n_p,N,n_snap);
snap.theta(:,:,1) = theta(:,:,1);
snap.k(:,:,:,:,1) = k;


%%% forward simulation
for i_t = 1:n-1
    for i_p = 1:n_p
        for i = 1:N
            dTheta = (omegaVect(i_p,i) + 1/N*sum(k(:,:,i_p,i).* sin(theta(:,:,i_t)-theta(i_p,i,i_t)),'all')) * dt; %+ zeta * sqrt(dt) * w(l,i);
            theta(i_p,i,i_t+1) = dTheta + theta(i_p,i,i_t);
        end
    end
    
    %%% updating indices of last spiking
    theta_t_w = wrapTo2Pi(squeeze(theta(:,:,i_t)));
    theta_tdt_w = wrapTo2Pi(squeeze(theta(:,:,i_t+1)));
    zeroXing_bool = theta_t_w - theta_tdt_w > minJump;
    last_spike_idx(zeroXing_bool) = i_t;
    
    %%% updating delta_t
    for i_p_pre = 1:n_p
        for i_pre = 1:N
            for i_p_post = 1:n_p
                for i_post = 1:N
                    if zeroXing_bool(i_p_pre,i_pre) || zeroXing_bool(i_p_post,i_post)
                        delta_t_val = (last_spike_idx(i_p_post,i_post) - ...
                            last_spike_idx(i_p_pre,i_pre)) * dt;
                        
                        %%% stdp rule
                        k_old = k(i_p_pre,i_pre,i_p_post,i_post);
                        if delta_t_val > 0
                            k_new = k_old + A_plus*exp(-delta_t_val/tau_plus);
                        elseif delta_t_val < 0
                            k_new = k_old - A_plus*beta*exp(-abs(delta_t_val)/(tau_R*tau_plus));
                        else
                            k_new = k_old;
                        end
                        
                        %%% k matrix update
                        k(i_p_pre,i_pre,i_p_post,i_post) = k_new;
                    end
                end
            end
        end
    end
    
    %%% updating time course of pop average couplings
    k_avg(:,:,i_t+1) = mean(k,[2 4]);%k(i_p_pre,i_pre,i_p_post,i_post);
    
    %%% updating snapshots if needed
    if mod(i_t,snapFact) == 0
        i_next = i_t/snapFact+1;
        snap.theta(:,:,i_next) = theta(:,:,i_t);
        snap.k(:,:,:,:,i_next) = k;
    end
    
end

%%% population order parameters
r_pop = sum(exp(1i*theta),2)/N;

end


