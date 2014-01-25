for i = 1:4
    figure(i)
    subplot(1,2,1);
    set(gca, 'YScale', 'linear');
    legend(gca, 'Location', 'NorthEast')
    subplot(1,2,2);
    set(gca, 'YScale', 'linear');
    legend(gca, 'Location', 'NorthEast')
    if i > 1
        export_fig -pdf -transparent -append
    else
        export_fig -pdf -transparent
    end
end
