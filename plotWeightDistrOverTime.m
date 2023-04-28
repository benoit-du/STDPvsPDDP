function [c_max_out,fname,counts,k_min_plot,k_max_plot] = ...
    plotWeightDistrOverTime(snap,k_avg,t,simPar,par,N,dirStr,titleStr,ft,doSave,c_max_in,k_min_plot,k_max_plot)

%%% 27-04-23    first revision
%%% Benoit Duchet, University of Oxford

% weights distribution over time
% snap.theta >> NaN(n_p,N,n_snap);
% snap.k >> NaN(n_p,N,n_p,N,n_snap);
n_k_edges = 100;
k_lim_pct = 5;
if isempty(k_min_plot)
k_min_plot = prctile(snap.k(:),k_lim_pct);
end
if isempty(k_max_plot)
k_max_plot = prctile(snap.k(:),100-k_lim_pct);
end
k_edges = linspace(k_min_plot,k_max_plot,n_k_edges);
yC = [mean([k_edges(2) k_edges(3)]) mean([k_edges(end-2) k_edges(end-1)])];
xC = [0 floor(simPar.n/simPar.snapFact)*simPar.snapFact*simPar.dt];
cax_cutoff = 95;

figure
for i_p_pre = 1:par.n_p
    for i_p_post = 1:par.n_p
        snaps = snap.k(i_p_pre,:,i_p_post,:,:);
        snaps = reshape(snaps,N*N,[]);
        counts = hist(snaps,k_edges);%% counts is n_edges x number of time points
        c_max = prctile(counts(:),cax_cutoff);
        
        subplot(par.n_p,par.n_p,sub2ind([par.n_p par.n_p],i_p_pre,i_p_post))
        colormap(jet)
        counts = counts(2:end-1,:);
        imagesc(xC,yC,counts)
        set(gca,'YDir','normal')
        c = colorbar;
        if isempty(c_max_in)
            caxis([0 c_max])
            c_max_out = c_max;
        else
            caxis([0 c_max_in])
            c_max_out = [];
        end
        ylabel(c,'counts','interpreter','latex');
        c.TickLabelInterpreter = 'latex';
        xlabel('time (s)','interpreter','latex')
        ylabel('$\kappa_{kl}$','interpreter','latex')
        
        hold on
        plot(t,squeeze(k_avg(i_p_pre,i_p_post,:)),'k','linewidth',1.5,'displayName',...
            ['N = ' num2str(par.N) ', avg pop ' num2str(i_p_pre) ' to ' num2str(i_p_post)])
        set(gca,'fontsize',ft)
    end
end

if doSave
    pause(1)
    fname = mySaveasFlex('dimXY',[12,11],'fNameNoNowStr',[dirStr filesep 'weightsEvol_' titleStr]);
else
    fname = [];
end
close


end

