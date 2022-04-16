% To toggle between PDF and CDF fit modes
function dataset = toggle_fit_mode(dataset,varargin)
    fprintf("Toggling fit mode for dataset %i", dataset.row)

    prefit_default = true;

    ip = inputParser;
    addRequired(ip,"dataset",@isstruct);
    addOptional(ip,"prefit",prefit_default,@islogical);
    parse(ip,dataset,varargin{:})

    prefit = ip.Results.prefit;

    x_all = dataset.bin_centres;
    C_all = dataset.CDF;
    P_all = dataset.PDF_normalised;
    [x,C,P] = get_valids(x_all,C_all,P_all);
    
    % Get prefit or random parameters for fit
    if prefit == true
        [M,S] = fit_lognormal_CDF_linear(x,C);
        p0 = [M,S];
    else
        p0 = rand(2,1);
    end

    % Do new fit, then calculate new values...
    if dataset.fit_by_CDF == true
        dataset.fit_by_CDF = false;
        [M,S] = fit_lognormal_PDF(x,P,p0);
    else
        dataset.fit_by_CDF = true;
        [M,S] = fit_lognormal_CDF(x,C,p0);
    end

    D43_fit = lognormal_moment(M,S,4,3);
    D32_fit = lognormal_moment(M,S,3,2);
    d10_fit = lognormal_quantile(M,S,10);
    d50_fit = lognormal_quantile(M,S,50);
    d90_fit = lognormal_quantile(M,S,90);
    phi = packing_fraction(M,S);
    mu = get_mu(M,S);
    sigma = get_sigma(M,S);
    COV = sigma/mu;

    % ...and replace values in dataset to be returned
    dataset.M = M;
    dataset.S = S;
    dataset.D43_fit = D43_fit;
    dataset.D32_fit = D32_fit;
    dataset.d10_fit = d10_fit;
    dataset.d50_fit = d50_fit;
    dataset.d90_fit = d90_fit;
    dataset.phi = phi;
    dataset.mu = mu;
    dataset.sigma = sigma;
    dataset.COV = COV;

end