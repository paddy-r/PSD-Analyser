## **Mathematical tools used by PSD Analyser**

### Hugh Patrick Rice, 2022

[![DOI](https://zenodo.org/badge/328645883.svg)](https://zenodo.org/badge/latestdoi/328645883)

## 1. Log-normal distribution

Many natural and manufactured particle species have size distributions that are log-normal, *i.e.* the logarithm of their size follows a normal distribution.

### Using the lognormal distribution in *PSD Analyser*

**1. Python**

**2. MATLAB**

**3. Microsoft Excel**

We note that the lognormal 

### Detailed description

## 2. Product difference algorithm (PDA)

For polydisperse particle species, *i.e.* those that have a distribution of sizes rather than a single one, there are many applications in which size fractions must be accounted for if a particle's physical properties depend on its size; for example, atmospheric and aerosol science and many other fluid-dynamic applications. The product difference algorithm is an example of a quadrature(-based) method of moments (QBMM, QMoM or QMM). It is an explicit, exact solution to the classical problem of moments and allows a discrete distribution to be expressed with an arbitrary number of elements with the same moments for computational efficiency.

### Using the PDA in *PSD Analyser*

**1. Python**

To obtain a lognormal fit to a size distribution, use either (1) *fit_lognormal_CDF* (takes sizes $d$ and CDF; $p_0$ is optional and is a set of lognormal parameters returned by *fit_lognormal_CDF_linear*), (2) *fit_lognormal_CDF_linear* ($d$, CDF) or (3) *fit_lognormal_PDF* ($d$, CDF, $p_0$ (optional)).

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
