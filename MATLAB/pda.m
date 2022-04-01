% HR 25/07/21
% To return N components of give distribution, where N is arbitrary
% using the product difference algorithm (PDA)

% HR 14/03/2017
% --
% To fit log-normal to particle size distribution
% Input should be row vector with d on first row, CDF on second row
% where d = particle size and CDF = cumulative dist. function, 0 < CDF < 1
% --
% CDF: C(d)  = 1/2 * [1 + erf((ln x - M)/S.sqrt(2))]
%      mu    = exp(M + S^2)
%      sigma = sqrt(exp(S^2 + 2M) * (exp(S^2) - 1))
% --
% Linearised:
%      ln d  = S.sqrt(2) * inv_err(2C - 1) + M
%      {__}                {_____________}
%       y    =     m              x        + c
% and S = m / sqrt(2)
% --
% PDF: P(d) = [1/(x.S.sqrt(2.pi)] * exp[-(ln x - M)^2/(2S^2)]
%
% --
% Quantile: Q(p) exp[M + sqrt(2.S^2) * inv_err(2p - 1)] 
%
% --
% HR 24/04/2017
% --
% To implement product difference algorithm (PDA) to model PSD as ternary
% (or N-ary) distribution (McGraw_71; Baker and Gammel_70 book;
% Marquisio_03; Mwasame_16 and Mwasame_16a; Gordon_68)
% --

clear all;

% All particle size data, pasted from Mastersizer data in Excel
% Data format is: size (micron); volume-weighted frequency; cumulative
% frequency
% d_h22     = [7.585776	8.709636	10	11.481536	13.182567	15.135612	17.378008	19.952623	22.908677	26.30268	30.199517	34.673685	39.810717	45.708819	52.480746	60.255959	69.183097	79.432823
% 0	0.007635	0.135005	0.310442	0.475434	0.558515	0.617157	0.910635	1.939477	4.28946	8.158915	12.980964	17.173874	18.739872	16.548617	11.479466	5.072049	0.602481
% 0	0.007635	0.14264	0.453082	0.928516	1.487031	2.104188	3.014823	4.9543	9.24376	17.402675	30.383639	47.557513	66.297385	82.846002	94.325468	99.397517	99.999998];

d_h22     = [8.709636	10	11.481536	13.182567	15.135612	17.378008	19.952623	22.908677	26.30268	30.199517	34.673685	39.810717	45.708819	52.480746	60.255959	69.183097	79.432823
0.007635	0.135005	0.310442	0.475434	0.558515	0.617157	0.910635	1.939477	4.28946	8.158915	12.980964	17.173874	18.739872	16.548617	11.479466	5.072049	0.602481
0.007635	0.14264	0.453082	0.928516	1.487031	2.104188	3.014823	4.9543	9.24376	17.402675	30.383639	47.557513	66.297385	82.846002	94.325468	99.397517	99.999998];

d_h16    = [34.673685	39.810717	45.708819	52.480746	60.255959	69.183097	79.432823	91.201084	104.712855	120.226443	138.038426	158.489319
0.012362	0.42421	2.158819	6.105918	12.005703	17.790159	20.449133	18.41713	12.907406	6.753966	2.47866	0.46569
0.012362	0.436572	2.595391	8.701309	20.707012	38.497171	58.946304	77.363434	90.27084	97.024806	99.503466	99.969156];

d_gb4060 = [2.187762	2.511886	2.884031	3.311311	3.801894	4.365158	5.011872	5.754399	6.606934	7.585776	8.709636	10	11.481536	13.182567	15.135612	17.378008	19.952623	22.908677	26.30268	30.199517	34.673685	39.810717	45.708819	52.480746	60.255959	69.183097	79.432823	91.201084	104.712855	120.226443	138.038426	158.489319	181.970086	208.929613	239.883292	275.42287	316.227766	363.078055	416.869383	478.630092	549.540874	630.957344	724.43596	831.763771	954.992586
0.011622	0.041607	0.061265	0.078719	0.095249	0.109825	0.125795	0.143752	0.16357	0.183986	0.203709	0.219662	0.229052	0.228285	0.214699	0.188279	0.149525	0.109847	0.060348	0.02189	0.022282	0.037277	0.060429	0.087608	0.108596	0.110986	0.089184	0.019803	0	0	0	0	0.054284	0.833449	2.326589	4.645038	7.519692	10.523236	12.936905	14.160764	13.80569	11.929445	8.994259	5.759525	2.804651
0.011622	0.053229	0.114494	0.193213	0.288462	0.398287	0.524082	0.667834	0.831404	1.01539	1.219099	1.438761	1.667813	1.896098	2.110797	2.299076	2.448601	2.558448	2.618796	2.640686	2.662968	2.700245	2.760674	2.848282	2.956878	3.067864	3.157048	3.176851	3.176851	3.176851	3.176851	3.176851	3.231135	4.064584	6.391173	11.036211	18.555903	29.079139	42.016044	56.176808	69.982498	81.911943	90.906202	96.665727	99.470378];

d_gb3040 = [316.227766	363.078055	416.869383	478.630092	549.540874	630.957344	724.43596	831.763771	954.992586	1096.478196	1258.925412	1445.439771	1659.586907
0.110684	1.056355	3.654673	8.166282	13.492893	17.56132	18.434032	15.786065	11.051514	6.344468	3.054888	1.131371	0.146073
0.110684	1.167039	4.821712	12.987994	26.480887	44.042207	62.476239	78.262304	89.313818	95.658286	98.713174	99.844545	99.990618];

d_calc   = [0.630957	0.724436	0.831764	0.954993	1.096478	1.258925	1.44544	1.659587	1.905461	2.187762	2.511886	2.884031	3.311311	3.801894	4.365158	5.011872	5.754399	6.606934	7.585776	8.709636	10	11.481536	13.182567	15.135612	17.378008	19.952623	22.908677	26.30268
0.10517	0.320736	0.654608	0.943462	1.296938	1.69634	2.17896	2.7447	3.391936	4.087646	4.799157	5.489189	6.119588	6.653002	7.048026	7.2733	7.295684	7.09573	6.674202	6.03487	5.23218	4.298399	3.333279	2.385141	1.536195	0.823786	0.255338	0.000183
0.10517	0.425906	1.080514	2.023976	3.320914	5.017254	7.196214	9.940914	13.33285	17.420496	22.219653	27.708842	33.82843	40.481432	47.529458	54.802758	62.098442	69.194172	75.868374	81.903244	87.135424	91.433823	94.767102	97.152243	98.688438	99.512224	99.767562	99.767745];

d_mh     = [0.954993	1.096478	1.258925	1.44544	1.659587	1.905461	2.187762	2.511886	2.884031	3.311311	3.801894	4.365158	5.011872	5.754399	6.606934	7.585776	8.709636	10	11.481536	13.182567	15.135612	17.378008	19.952623	22.908677	26.30268	30.199517	34.673685	39.810717	45.708819	52.480746	60.255959	69.183097	79.432823
0.240277	1.179956	2.11795	3.265372	4.430556	5.536258	6.461619	7.14195	7.551812	7.698655	7.59965	7.278433	6.769304	6.132137	5.410137	4.664515	3.916301	3.213404	2.558128	1.984723	1.483509	1.071384	0.738507	0.488814	0.3107	0.195809	0.128745	0.095176	0.081095	0.076037	0.068625	0.054399	0.033771
0.240277	1.420233	3.538183	6.803555	11.234111	16.770369	23.231988	30.373938	37.92575	45.624405	53.224055	60.502488	67.271792	73.403929	78.814066	83.478581	87.394882	90.608286	93.166414	95.151137	96.634646	97.70603	98.444537	98.933351	99.244051	99.43986	99.568605	99.663781	99.744876	99.820913	99.889538	99.943937	99.977708];

d_bar    = [0.630957	0.724436	0.831764	0.954993	1.096478	1.258925	1.44544	1.659587	1.905461	2.187762	2.511886	2.884031	3.311311	3.801894	4.365158	5.011872	5.754399	6.606934	7.585776	8.709636	10	11.481536	13.182567	15.135612	17.378008	19.952623	22.908677	26.30268	30.199517	34.673685	39.810717	45.708819	52.480746	60.255959	69.183097	79.432823	91.201084	104.712855	120.226443	138.038426	158.489319	181.970086	208.929613	239.883292	275.42287	316.227766	363.078055	416.869383	478.630092	549.540874	630.957344	724.43596	831.763771	954.992586
0.101416	0.235209	0.387333	0.471715	0.548949	0.618211	0.710502	0.844474	1.035343	1.287211	1.607949	2.009277	2.500496	3.085807	3.75555	4.502314	5.283456	6.060006	6.751724	7.300573	7.616659	7.646689	7.35361	6.730918	5.840276	4.750634	3.600736	2.492736	1.551926	0.815274	0.31155	0.092843	0.048566	0.044578	0.051166	0.057903	0.059103	0.053749	0.043289	0.029314	0.010179	0	0	0	0	0	0	0.021801	0.08241	0.150301	0.219117	0.280385	0.322231	0.32987
0.101416	0.336625	0.723958	1.195673	1.744622	2.362833	3.073335	3.917809	4.953152	6.240363	7.848312	9.857589	12.358085	15.443892	19.199442	23.701756	28.985212	35.045218	41.796942	49.097515	56.714174	64.360863	71.714473	78.445391	84.285667	89.036301	92.637037	95.129773	96.681699	97.496973	97.808523	97.901366	97.949932	97.99451	98.045676	98.103579	98.162682	98.216431	98.25972	98.289034	98.299213	98.299213	98.299213	98.299213	98.299213	98.299213	98.299213	98.321014	98.403424	98.553725	98.772842	99.053227	99.375458	99.705328];

d_sinclair_62_f2 = [147	125	104	88	74
31.7	49.8	12.2	0.3	0.1
94.1	62.4	12.6	0.4	0.1];

d_smith_55_f9    = [420	295	210	150	104
0.9	28	64.5	6	0.5
99.999	99	71	6.5	0.5];

% Choose distribution; convert from percentages
d      = d_h22;
bins = [0 d(1,:)];
PDF = d(2,:)/100;
CDF = d(3,:)/100;

PDF_sum = sum(PDF);

% Recalculate bin-median particle sizes; truncate first column (as
% Mastersizer gives frequency dist "b/t" sizes)
for i = 1:(length(bins)-1)
    x(i) = mean(bins(i:(i+1)));
    bin_widths(i) = bins(i+1) - bins(i);
    PDF_corr(i) = PDF(i)/(bin_widths(i) * PDF_sum);
%     CDF_corr(i) = CDF(i)/(bin_widths(i) * PDF_sum)
end

% d(5,:) = [0 diff(d(1,:))];
% d      = d(:,2:length(d));

% %%
% % Find jamming fraction (i.e. high-shear packing fraction)
% phi_rcp = 0.644;
% A       = (m(2)*m(4)/m(3)^2) * phi_rcp/(1-phi_rcp);
% phi_j   = A/(1+A);

% %% Ternary viscosity (Mwasame_16, 16a)
% phi_total = 0.1;
% phi       = w.*phi_total;

% % Find alpha_star, beta_1 and beta_2 for Mwasame_16a model
% K          = 100; % Given as limit to infinity in Mwasame_16a
% alpha_star = D(2) + (D(3)-D(2))*((phi(1) + phi(3))/(phi(1) + phi(2) + phi(3)))^K;
% fprintf('K = %d\nalpha_star = %d\n',K,alpha_star);

% % Find beta_1 and beta_2 by iteration
% beta_2 = 0.9;
% for i = 1:10
%     beta_1     = beta(D(1), alpha_star, phi(1), phi(2) + beta_2*phi(3));
%     beta_2     = beta(alpha_star, D(3), beta_1*phi(1) + phi(2), phi(3));
%     fprintf('beta_1 = %d\nbeta_2 = %d\n',beta_1,beta_2);
% end

%% Log-normal and normal fit calculation and plots

% Log-normal fit to volume data, linearised
X = erfinv(2*CDF - 1);
Y = log(x);

% Linear fit
p_logn     = polyfit(X,Y,1);
S     = p_logn(1)/sqrt(2);
M     = p_logn(2);
mu_logn    = exp(M + S^2/2);
sigma_logn = sqrt(exp(S^2 + 2*M) * (exp(S^2) - 1));

% HR 20/03/22 Fit via optmisation
B = fminsearch(@(b) norm(CDF - logncdf(x,b(1),b(2))), rand(2,1));
M = B(1);
S = B(2);
PDF_normalised = PDF./(bin_widths*sum(PDF));

% HR 26/06/2018
PDF_fit   = (1./(S.*sqrt(2*pi)*x(:))).*exp(-((log(x(:))-M).^2)/(2*(S^2)));
CDF_fit   = 0.5*(1+erf((log(x(:))-M)/(S*sqrt(2))));

% Goodness of fit, log-normal case
y_fit  = polyval(p_logn,X);
y_res  = Y - y_fit;
ss_res = sum(y_res.^2);
ss_tot = (length(Y)-1) * var(Y);
r_sq_logn   = 1 - ss_res/ss_tot;

% Normal fit to volume data, linearised
p_norm   = polyfit(X,x,1);
sigma_norm = p_norm(1)/sqrt(2);
mu_norm    = p_norm(2);

% Goodness of fit, log-normal case
y_fit   = polyval(p_norm,X);
y_res   = x - y_fit;
ss_res  = sum(y_res.^2);
ss_tot  = (length(x(:))-1) * var(x(:));
r_sq_norm = 1 - ss_res/ss_tot;

sprintf('R-squared for log-normal case: %f',r_sq_logn);
sprintf('R-squared for normal case:     %f',r_sq_norm);

N = 3;
[D,w,m] = PDA(x,PDF,N);
% D = PDA_results(1,:);
% w = PDA_results(2,:);

for i = 0:(2*N-1)
    m_fit(i+1) = sum(D.^i.*w);
    m_ln(i+1)  = exp(i*M + 0.5*i^2*S^2);
end

% Normalise for comparison
for i = 1:length(m)
    m_norm(i)      = m(i)/m(2)^i;
    m_fit_norm(i)  = m_fit(i)/m(2)^i;
    m_ln_norm(i)   = m_ln(i)/m(2)^i;
    rel_error(i)   = abs((m(i)-m_fit(i))/m(i));
    rel_error2(i)  = abs((m(i)-m_ln(i))/m(i));
    fprintf('m(%i) = %d, m_fit(%i) = %d, m_ln(%i) = %d\n',i-1,m(i),i-1,m_fit(i),i-1,m_ln(i));
    fprintf('Rel. error vs. PDA = %d, vs. log-normal fit = %d\n',rel_error(i),rel_error2(i));
end

%%
sz   = 12;
font = 'Cambria';

% (1) Plot lognormal fit in linearised form
hold off
plot(X,Y,'blackx')
hold all
a = xlim;
b = ylim;
c = linspace(a(1),a(2),2);
plot(c,(c.*p_logn(1,1) + p_logn(1,2)),'black--');
hold off

axis([a(1) a(2) b(1) b(2)]);
legend('Sizing data','Linear fit','Location', 'SE');
xlabel('{\itS}\surd2 erf^{-1}(2{\itC}-1) + \itM','FontSize',sz,'FontName',font);
ylabel('ln \itd','FontSize',sz,'FontName',font);
set(gca,'FontSize',sz,'FontName',font);

pause;

% (2) Plot normal fit in linearised form
hold off

plot(X,x,'blackx')
hold all

a = xlim;
b = ylim;
c = linspace(a(1),a(2),2);
plot(c,(c.*p_norm(1) + p_norm(2)),'black--');
hold off

axis([a(1) a(2) b(1) b(2)]);
legend('Sizing data','Linear fit','Location', 'SE');
xlabel('\sigma\surd2 erf^{-1}(2{\itC}-1) + \it\mu','FontSize',sz,'FontName',font);
ylabel('\itd','FontSize',sz,'FontName',font);
set(gca,'FontSize',sz,'FontName',font);

pause;

% Get plot width
d_min   = min(x);
d_max   = max(x);
d_range = d_max - d_min;
zoom    = 0.1;
d_min   = d_min - zoom*d_range;
d_max   = d_max + zoom*d_range;

% (3) Plot CDF
plot(x,CDF,'blackx');
hold all

x_plot = linspace(d_min,d_max);
plot(x_plot,logncdf(x_plot,M,S),'black--');
plot(x_plot,normcdf(x_plot,mu_norm,sigma_norm),'black');

plot([d_min d_max],[0.1 0.1],'black:');
plot([d_min d_max],[0.5 0.5],'black:');
plot([d_min d_max],[0.9 0.9],'black:');
hold off

axis([d_min d_max 0 1]);
legend('Sizing data','Log-normal fit','Normal fit','10th, 50th, 90th %ile','Location', 'Best');
xlabel('Particle size (\mum)','FontSize',sz,'FontName',font);
ylabel('CDF','FontSize',sz,'FontName',font);
set(gca,'FontSize',sz,'FontName',font);

pause;

% (4) Plot PDF vs. d
plot(x,PDF_normalised,'blackx');
hold all

x_plot = linspace(d_min,d_max);
y_pdf = lognpdf(x_plot,M,S);
plot(x_plot,y_pdf,'black--');
% plot(x,PDF_fit,'black--');
plot(x_plot,normpdf(x_plot,mu_norm,sigma_norm),'black');

% Get d values corresp. to 10th, 50th and 90th percentiles
d_10 = quantile(M,S,0.1);
d_50 = quantile(M,S,0.5);
d_90 = quantile(M,S,0.9);

y_max = max(y_pdf)*1.1;
axis([d_min d_max 0 y_max]);

b = ylim;
plot([d_10 d_10],[b(1) b(2)],'black:');
plot([d_50 d_50],[b(1) b(2)],'black:');
plot([d_90 d_90],[b(1) b(2)],'black:');
hold off

legend('Sizing data','Log-normal fit','Normal fit','10th, 50th, 90th %ile','Location', 'Best');
xlabel('Particle size (\mum)','FontSize',sz,'FontName',font);
ylabel('PDF','FontSize',sz,'FontName',font);
set(gca,'FontSize',sz,'FontName',font);

x = [1 2 3 4 5];
P = [0.1 0.2 0.3 0.25 0.15];
[x_,P_,m_] = PDA(x,P,3);

% HR 12/03/22 PDA function
function [D,w,m] = PDA(x,PDF,N)
    N = 3;
    for k = 1:(2*N)
        m(k) = sum(PDF(:).*x(:).^(k-1));
    end
    
    P = zeros(2*N+1);
    P(1,1) = 1;
    for i = 1:(2*N)
        P(i,2) = (-1)^(i-1) * m(i);
    end
    
    for j = 3:(2*N+1)
        for i = 1:(2*N+2-j)
            P(i,j) = P(1,j-1)*P(i+1,j-2) - P(1,j-2)*P(i+1,j-1);
        end
    end

    alpha(1) = 0;
    for i = 2:2*N
        alpha(i) = P(1,i+1)/(P(1,i)*P(1,i-1));
    end
    
    for i = 1:N
        a(i) = alpha(2*i) + alpha(2*i - 1);
    end
    
    for i = 1:(N-1)
        b(i) = -sqrt(alpha(2*i + 1)*alpha(2*i));
    end
    
    J = gallery('tridiag',b,a,b)
    [w,D] = eig(full(J))
    
    D = diag(D);
%     w = w(1,:).^2;
%     w = w';
    w = w(:,1).^2

end

% Get x corresp. to pth quantile of lognormal distribution
function x_p = quantile(M,S,p)
    x_p = exp(M + sqrt(2*S^2)*erfinv(2*p-1.0));
end