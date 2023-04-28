function useLatex(bool)
%%% 06-06-18    Switching between latex and default interpreter for figures
%%%             Has to be called before plotting the figure
%%%             Then got to use $ for math (not needed with tex)

if bool
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
        set(groot, 'defaultLegendInterpreter','latex');
else
        set(groot, 'defaultAxesTickLabelInterpreter','tex'); 
        set(groot, 'defaultLegendInterpreter','tex');
end

end
    