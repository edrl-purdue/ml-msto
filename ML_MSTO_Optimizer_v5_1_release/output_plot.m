% OUTPUT ITERATION INFO AND PLOTS
function output_plot(xMacro,xWeight,macro,xMicro,c,loop,change)

    %% PRINT RESULTS
    fprintf('     |     Obj.:%.4e Vol.:%7.3f ch.:%7.3f      |\n\n',c,mean(xMacro(:)),change);

    %% PLOT RESULTS
    if macro.displayflag
        xMacro = reshape(xMacro,[macro.nelx, macro.nely, max(macro.nelz,1)]);
        figure(1)
        clf;
        display_top(real(xMacro)); drawnow;
        
        if macro.alg ~= 1
            figure(2)
            h = findobj(gca,'Type','line');
            if ~isempty(h)
                plot([h.XData, loop], [h.YData, c],'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','m','Marker','diamond'); xlabel('Iterations'); ylabel('Function value'); drawnow;
            else
                plot(loop, c,'LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','m','Marker','diamond'); xlabel('Iterations'); ylabel('Function value'); drawnow;
            end
        end

        if all(macro.phys == [1 2])
            xWeight = reshape(xWeight,[macro.nelx, macro.nely, max(macro.nelz,1)]);
            figure(4)
            clf;
            display_obj(xWeight); drawnow;
        end

        if ~isempty(cell2mat(xMicro))
            figure(3)
            clf;
            display_top(cell2mat(reshape(xMicro,[macro.nelx, macro.nely, max(macro.nelz,1)]))); drawnow;
        end
    end
end