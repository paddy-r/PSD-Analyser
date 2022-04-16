% Plot (normalised) PDF, optionally with fit curve
function plot_PDF(x,P,varargin)
    fit_data_default = [];
    log_mode_default = false;
    x_label_default = "Particle diameter (\mum)";
    y_label_default = "PDF (normalised)";

    ip = inputParser;
    addRequired(ip,'x',@(x) numel(x)>0);
    addRequired(ip,'P',@(P) numel(P)>0);
    addOptional(ip,'fit_data',fit_data_default);
    addOptional(ip,'log_mode',log_mode_default);
    addOptional(ip,'x_label',x_label_default,@isstring);
    addOptional(ip,'y_label',y_label_default,@isstring);
    parse(ip,x,P,varargin{:});

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
        semilogx(x,P,'blackx');
    else
        plot(x,P,'blackx');
    end
    
    if numel(fit_data) > 0
        hold all
        M = fit_data(1)
        S = fit_data(2)
        if log_mode == true
            x_plot = linspace(min(x),max(x),1000);
            y_pdf = lognormal_PDF(x_plot,M,S);
            semilogx(x_plot,y_pdf,'black--');
        else
            x_plot = linspace(min(x),max(x),100);
            y_pdf = lognormal_PDF(x_plot,M,S);
            plot(x_plot,y_pdf,'black--');
        end
    end
    
    if x_label ~= ""
        xlabel(x_label);
    end
    if y_label ~= ""
        ylabel(y_label);
    end
    hold off
    
end