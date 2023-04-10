function plotfn(fn, xlim, ylim)
    % plotting function
    figure;
    [ptsx, ptsy] = meshgrid(linspace(xlim(1),xlim(2),800), linspace(ylim(1),ylim(2),800));
    fn_vals = fn(ptsx,ptsy);
    [min_val] = min(fn_vals,[], "all");
    [max_val] = max(fn_vals,[], "all");
    mesh(ptsx ,ptsy, fn_vals);
    colormap default; colorbar
    hold on
    plot3(ptsx(fn_vals == min_val),ptsy(fn_vals == min_val), fn_vals(fn_vals == min_val), "xr");
    plot3(ptsx(fn_vals == max_val),ptsy(fn_vals == max_val), fn_vals(fn_vals == max_val), "xg");
    title("f(x)");
    xlabel("x","FontWeight", "bold"); ylabel("y", "FontWeight", "bold"); zlabel("f(x)", "FontWeight", "bold");
    legend("function", "Min value", "Max value", "Location", "best");
end