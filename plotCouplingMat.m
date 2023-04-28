function [fname,kmat] = plotCouplingMat(snap,N,n_t_plot,sortType,omegaVect,idx_phVel_sort,dirStr,titleStr,ft,clims)

%%% 27-04-23    first revision
%%% Benoit Duchet, University of Oxford

%%% sortType: 0 is no sorting, 1 by natural freq, 2 by phase velocity

% coupling matrix over time
% snap.k >> NaN(n_p,N,n_p,N,n_snap);
% k(i_p_pre,i_pre,i_p_post,i_post);
n_t = size(snap.theta,3);
i_t_plot = round(linspace(1,n_t,n_t_plot));

figure
for i_t_ind = 1:4
    i_t = i_t_plot(i_t_ind);
   
    switch sortType
        case 0
            %%% no sorting
            kmat = squeeze(snap.k(1,:,1,:,i_t));
        case 1
            %%% sorting by natural frequency
            [~,idx_freq_sort] = sort(omegaVect(1,:));
            kmat = squeeze(snap.k(1,idx_freq_sort,1,idx_freq_sort,i_t));
            
        case 2
            %%% sorting by phase velocity at steady state
            kmat = squeeze(snap.k(1,idx_phVel_sort,1,idx_phVel_sort,i_t));
    end
    
    imagesc([1 N],[1 N],kmat)
    colormap(jet)
    set(gca,'YDir','normal')
    c = colorbar;
    ylabel(c,'$\kappa_{kl}$','interpreter','latex');
    c.TickLabelInterpreter = 'latex';
    xlabel('post index (k)','interpreter','latex')
    ylabel('pre index (l)','interpreter','latex')
    set(gca,'fontsize',ft)
    if ~isempty(clims)
        caxis(clims)
    end
    pause(1)
    fname = mySaveasFlex('dimXY',[12,11],'fNameNoNowStr',[dirStr filesep titleStr '_t' num2str(i_t_ind)]);
    close
end


end
