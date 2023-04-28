function f = getFourierEqRule(eps,A_p,A_d,tau_p,tau_d,n,verbose)

%%% 27-04-23    first revision
%%% Benoit Duchet, University of Oxford

n_pts = 500;
phi = linspace(0,2*pi,n_pts);
dK = eps*(A_p*exp(-phi/tau_p)-A_d*exp((phi-2*pi)/tau_d));

ftStr = 'a0';
for i = 1:n
    ftStr = [ftStr ' + a' num2str(i) '*cos(' num2str(i) '*x) + b' num2str(i) '*sin(' num2str(i) '*x)'];
end
ft = fittype(ftStr);
f = fit(phi(:), dK(:), ft, 'Algorithm', 'Levenberg-Marquardt');

if verbose
    figure
    subplot(1,2,1)
    plot(f,phi,dK)
    xlabel('$\varphi$','interpreter','latex')
    ylabel('$\dot{K}_{ij}$','interpreter','latex')
    title([num2str(n) ' Fourier modes'])
    
    subplot(1,2,2)
    phi_broad = linspace(-pi,3*pi,2*n_pts);
    plot(phi_broad,feval(f,phi_broad))
    xlabel('$\varphi$','interpreter','latex')
    title('showing periodicity')
end

end