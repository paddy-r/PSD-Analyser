% Plot CDF, optionally with fit curve
% Args (required):
%   x           particle diameter (or x value)
%   C           cumulative distribution function, PDF
% Args (optional):
%   fit_data    M,S (lognormal parameters) for plotting fit curve
%   log_mode    if false (default) print log-x scale; linear if true
%   x_label     label text for x axis
%   y_label     label text for y axis
function current_figure = plot_PDF(x,C,varargin)
    fit_data_default = [];
    log_mode_default = false;
    x_label_default = "Particle diameter (\mum)";
    y_label_default = "CDF";

    ip = inputParser;
    addRequired(ip,'x',@(x) numel(x)>0);
    addRequired(ip,'C',@(C) numel(C)>0);
    addOptional(ip,'fit_data',fit_data_default);
    addOptional(ip,'log_mode',log_mode_default);
    addOptional(ip,'x_label',x_label_default,@isstring);
    addOptional(ip,'y_label',y_label_default,@isstring);
    parse(ip,x,C,varargin{:});

    sprintf("%s",string(ip.Results.fit_data));
    sprintf("%s",string(ip.Results.log_mode));
    sprintf("%s",ip.Results.x_label);
    sprintf("%s",ip.Results.y_label);

    log_mode = ip.Results.log_mode
    fit_data = ip.Results.fit_data;
    x_label = ip.Results.x_label;
    y_label = ip.Results.y_label;

    % Plot PSD data
    hold off
    if log_mode == true
        semilogx(x,C,'blackx');
    else
        plot(x,C,'blackx');
    end
    
    if numel(fit_data) > 0
        hold all
        M = fit_data(1)
        S = fit_data(2)
        if log_mode == true
            x_plot = logspace(log10(min(x)),log10(max(x)),100);
            y_cdf = lognormal_CDF(x_plot,M,S);
            semilogx(x_plot,y_cdf,'black--');
        else
            x_plot = linspace(min(x),max(x),100);
            y_cdf = lognormal_CDF(x_plot,M,S);
            plot(x_plot,y_cdf,'black--');
        end
    end
    
    if x_label ~= ""
        xlabel(x_label);
    end
    if y_label ~= ""
        ylabel(y_label);
    end
    hold off
    
    % Grab figure handle to be returned
    current_figure = gcf
end