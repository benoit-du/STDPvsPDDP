%%% 27-04-23    first revision
%%% Benoit Duchet, University of Oxford

clearvars
close all

timeInS         = 150;           %%% duration of simulation (biological time in seconds)
dt              = 1E-3;             %%% simulation time step (in seconds)

N               = 60;               %%% number of oscillators
Omega           = 2*pi*4.96;        %%% average natural frequency of oscillators (rd/s)
Delta           = 0.6*pi;           %%% standard deviation of oscillator frequency distribution (rd/s)

tau_plus        = 0.0168;           %%% timescale of long-term potentiation (LTP) in seconds
tau_R           = 2;                %%% = tau_minus/tau_plus, i.e. ratio of the timescale of long-term depression (LTD) to the timescale of LTP 
A_plus          = 0.2;              %%% magnitude of LTP
beta            = 0.5;              %%% = A_minus/A_plus, i.e. ratio of the magnitude of LTD to LTP    
n_modes         = 5;                %%% number of Fourier modes in the PDDP approximation

theta_0_std     = pi/3;             %%% standard deviation of the distribution of initial phases 
kappa_0         = 12;               %%% mean of the distribution of initial coupling strengths 
kappa_0_std     = 0.2;              %%% standard deviation of the distribution of initial coupling strength

outputDir       = 'Fig';            %%% directory where output will be saved
seed            = floor(rand*1E6);  %%% seed used to initialise the random number generator (set to a fixed value for repeatable results)

tic
mkdir(outputDir)
STDPvsPDDP(N,Omega,Delta,theta_0_std,kappa_0,kappa_0_std,tau_plus,tau_R,beta,A_plus,n_modes,outputDir,seed,timeInS,dt);
STDPvsEbPDDP(N,Omega,Delta,theta_0_std,kappa_0,kappa_0_std,tau_plus,tau_R,beta,A_plus,n_modes,outputDir,seed,timeInS,dt);
toc
