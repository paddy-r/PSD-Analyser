## **Mathematical tools used by PSD Analyser**

### Hugh Patrick Rice, 2022

[![DOI](https://zenodo.org/badge/328645883.svg)](https://zenodo.org/badge/latestdoi/328645883)

## 1. Log-normal distribution

Many natural and manufactured particle species have size distributions that are log-normal, *i.e.* the logarithm of their size follows a normal distribution.

### Using the lognormal distribution in *PSD Analyser*

**1. Python**

To obtain a lognormal fit to a size distribution, use either (1) *fit_lognormal_CDF* (takes sizes $d$ and CDF; $p_0$ is optional and is a set of lognormal parameters returned by *fit_lognormal_CDF_linear*), (2) *fit_lognormal_CDF_linear* ( $d$, CDF) or (3) *fit_lognormal_PDF* ( $d$, CDF, $p_0$ (optional)).

To obtain the value of the PDF or CDF of a lognormal distribution with parameters $M$ and $S$ at a size $x$, use (1) *lognormal_PDF* or (2) *lognormal_CDF*, respectively.

**2. MATLAB**

The MATLAB functions are named exactly the same as the Python ones, but differ in that $p_0$ is a required argument in *fit_lognormal_CDF* and *fit_lognormal_PDF*.

**3. Microsoft Excel**

We note that the lognormal function is not implemented in Microsoft Excel, but that the related inverse error function, $\Gamma^{–1}$ (see Equation (3)), is. Details are given at the end of the next section.

### Detailed description

The lognormal cumulative distribution function (CDF), $C$, is:

$$ C(d) = \frac{1}{2} \left[ 1+\mathrm{erf} \left( \frac{\ln⁡(d-M)}{S\sqrt2} \right) \right], \tag{1} $$

where $d$ is the particle size in microns and $M$ and $S$ are the log-normal parameters and $\mathrm{erf}$ is the error function. Equation (1) can be linearised as follows in order to fit measured data:

$$ \mathrm{erf}^{-1} \left[ 2C(d)-1 \right] = \frac{\ln⁡d}{S\sqrt2} + \frac{M}{S\sqrt2} \tag{2}, $$

*i.e.* $S$ and $M$ can be found from the gradient and intercept, respectively, of a plot of the LHS of Equation (2) *vs.* $\ln d$. The $q$ th quantile is given by:

$$ d_q = \exp \left[ M + \sqrt{2S^2} \mathrm{erf}^{-1} \left( 2q-1 \right) \right], \tag{3} $$

where $\mathrm{erf}^{–1}” is the inverse error function, and the median value – *i.e.* $q = 50$ or $d = d_{50}$ – is $\exp(M)$. The $n$ th moment of a log-normal distribution is:

$$ m_n = \exp \left( nM + \frac{1}{2} n^2 S^2 \right), \tag{4} $$

and particle size metrics are commonly given in the form $d[a,b] = m_a/m_b$, for example the volume-weighted mean, $d[4,3]$, and the surface-weighted or Sauter mean, $d[3,2]$.

There is no built-in inverse error function in Microsoft Excel, but the inverse gamma function, $\Gamma^{–1}$ (implemented as ```GAMMAINV(𝑥, α, β)``` or newer version ```GAMMA.INV(𝑥, α, β))```, is provided and is related to the inverse error function as follows, with $\alpha = 0.5$ and $\beta = 1$, i.e.:

$$ \sqrt{\Gamma^{-1} (x,0.5,1)} = \mathrm{erf}^{-1} (x). \tag{5} $$

The first argument to Excel’s inverse gamma function must be in the range $0 < x < 1$, which means that only the lower half of the particle size data can be used for the fit. However, since the inverse error function is odd – *i.e.* $\mathrm{erf}^{–1}(–x) = –\mathrm{erf}^{–1}(x)$ – the full particle size distribution can be accessed – *i.e.* $–1 < x < 1$ – using a conditional statement in Excel as follows:

```
=IF(cell>0,SQRT(GAMMAINV(cell,0.5,1)),-SQRT(GAMMAINV(–cell,0.5,1))) 
```

where cell indicates a reference to the relevant cell in the spreadsheet (in this case, the cell containing $x = 2C – 1$). A similar conditional formula can be written to evaluate quantiles *via* Equation (3).

## 2. Product difference algorithm (PDA)

For polydisperse particle species, *i.e.* those that have a distribution of sizes rather than a single one, there are many applications in which size fractions must be accounted for if a particle's physical properties depend on its size; for example, atmospheric and aerosol science and many other fluid-dynamic applications. The product difference algorithm is an example of a quadrature(-based) method of moments (QBMM, QMoM or QMM). It is an explicit, exact solution to the classical problem of moments and allows a discrete distribution to be expressed with an arbitrary number of elements with the same moments for computational efficiency.

### Using the PDA in *PSD Analyser*

**1. Python**

**2. MATLAB**

### Detailed description

The $k$ th statistical moment, $m_k$, of a distribution of a variable $x$ is:

$$ m_k = \int x^k f(x)dx \approx \sum_{i=1}^n w_i x_i^k, \ k = 0,\dots,N, \tag{1} $$
	

where $f(x)$ is the distribution function in the range $x$ to $x$ + $dx$, and $w_i$ and $x_i$ are the $i$th weights and abscissas of a corresponding discrete distribution; $w_i$ may alternatively be expressed as $W_i$ $\Delta x_i$ in the case of a measured distribution with variable bin sizes, where $W_i$ is a suitably modified weight. For closure, it is required that $N = 2n-1$, and $m_0$ = 1 if the distribution is suitably normalised.

The question of whether a function $f$ (or $w_i$ and $x_i$ in the discrete case) exists that satisfies Equation (1), where the values of $m_k$ are known, is a form of the classical problem of moments (Shohat and Tamarkin, 1943; Wheeler and Gordon, 1971). Equation (1) or similar expressions occur widely in modelling of particulate kinetics (Marchisio *et al.*, 2003; McGraw, 1997) and modelling of a distribution (particle size, for example) as a small number of moments offers computational advantages. In the discrete case Equation (1) produces $N+1$ non-linear equations the solution of which is not straightforward. An explicit, $O(n^2)$ method for solving Equation (1) for $w_i$ and $x_i$ for arbitrary values of $N$ was described by Gordon (1968) and is referred to as the product-difference algorithm to differentiate it from other, similar ones such as quotient-difference algorithms (Rutishauser, 1957). We are not aware of a complete, correct description of the PDA in recipe form, including indices, for numerical computation and what follows is presented as such.

The procedure for implementation the PDA begins with the construction of a $2N+1$ square matrix $P$ such that $P_{1,1} = 1$ and $P_{i,1} = 0$ otherwise, and:

$$ P_{i,2} = (-1)^{i-1} m_{i-1}, \ i = 1,\dots,2n, \tag{2} $$
	
$$ P_{i,j} = P_{1,j-1} P_{i+1,j-2}-P_{1,j-2} P_{i+1,j-1}, \ i = 1,\dots,2n+2-j, \ j = 3,\dots,2n+1. \tag{3} $$
	
Explicit expressions and tabulations of $P$ for various values of $n$ and $m_k$ (Gordon, 1968; Marchisio *et al.*, 2003; Wheeler and Gordon, 1971). A $2n$-sized vector, $\alpha$, is then constructed such that $\alpha_1 = 0$ and:

$$ α_i=P_{1,i+1}/(P_{1,i} P_{1,i-1} ), \ i = 2,\dots,2n. \tag{4} $$	

Finally, $n$- and $(n-1)$-sized vectors, $a$ and $b$ respectively, are constructed such that:

$$ a_i = α_2i + α_{2i-1}, \ i = 1,\dots,n, \tag{5} $$
	
$$ b_i = -\sqrt{α_{2i+1} + α_2i}, \ i = 1,\dots,n-1, \tag{6} $$

and are then used to construct a tridiagonal matrix, $J$, with a forming the main diagonal and $b$ both the sub- and superdiagonals. The eigenvalues of $J$ give the abscissas, $x_i$, directly and the first components, $v_{i,1}$, of each eigenvector of $J$ give the weights, $w_i$, as follows (McGraw, 1997):

$$ w_i=m_0 v_{i,1}^2, \ i = 1,\dots,n. \tag{7} $$

The PDA was used recently and very successfully by Mwasame *et al.* (2016) to model volume-weighted particle size distributions (with $n = 3$) in order to implement a model of the viscosity of multiphase (solid-liquid) mixtures with arbitrary solid-phase size distributions. In the case of $n = 3$, *i.e.* a ternary distribution, Equation (1) becomes:

$$ m_k = \sum_{i=1}^3 w_i x_i^k = w_1 x_1^k + w_2 x_2^k + w_3 x_3^k, \ k=0,\dots,5, \tag{8} $$

and $a$ and $b$ (Equations (5) and (6)) form the matrix $J$ as follows (McGraw, 1997):

$$ J = \begin{pmatrix}
  b_1 & a_1 & \\
  a_1 & b_2 & a_2 \\
      & a_2 & b_3 \\
      \end{pmatrix}. \tag{9} $$
