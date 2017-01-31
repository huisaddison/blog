Title: Discussion: Optimal Rates of Convergence for Covariance Matrix Estimation, Part 2
Date: 2017-01-30
Modified: 2017-01-30
Category: Statistics
Tags: math, stats, covariance, minimax, risk
Slug: discuss-covmat-optimal-convergence-pt2
Summary: Discussion of paper on minimax estimation of covariance matrices, continued.

In this post, I will continue my discussion of [Optimal Rates of Convergence
for Covariance Matrix Estimation](https://arxiv.org/abs/1010.3866).  A
discussion of their first result (showing an upper bound on the risk of their
proposed tapering estimator) can be found in a [previous
post]({filename}./discuss-covmat-optimal-convergence-pt1.md).

## Minimax Lower Bound
We will review their proof of a minimax lower bound among all estimators that
matches the upper bound previously proved for their tapering estimator.

**Theorem 3.**  _Suppose $p \leq \exp(\gamma n)$ for some constant
$\gamma > 0$.  The minimax risk for estimating the covariance matrix
$\Sigma$ over $\mathcal{P}_\alpha$ under the operator norm satisfies_:
\begin{align*}
  \inf_{\hat\Sigma}\sup_{\mathcal{P}_\alpha}
  \mathbf{E}\Vert\hat\Sigma - \Sigma\Vert^2
  \geq
  cn^{-\frac{2\alpha}{2\alpha + 1}} + c\frac{\log p}{n} 
\end{align*}

In the words of the authors,
> _The basic strategy underlying the proof of Theorem 3 is to carefully
  construct a finite collection of multivariate normal distributions and
  calculate the total variation affinity between pairs of probability measures
  in the collection._

## Parameter Space Specification
Oftentimes, the proof of a minimax lower bound is accompanied with the
specification of a smaller parameter space that is easier to analyze.
Intuitively, the restriction of the parameters to a subspace is permissible as
for $\mathcal{F}'\subset \mathcal{F}$:
\begin{align*}
    \sup_{\Sigma\in\mathcal{F}'} R(\Sigma, \hat\Sigma)
    \leq
    \sup_{\Sigma\in\mathcal{F}} R(\Sigma, \hat\Sigma)
\end{align*}
for all $\hat\Sigma$.  It follows, then, that:
\begin{align*}
    \inf_{\hat\Sigma}\sup_{\Sigma\in\mathcal{F}'} R(\Sigma, \hat\Sigma)
    \leq
    \inf_{\hat\Sigma}\sup_{\Sigma\in\mathcal{F}} R(\Sigma, \hat\Sigma)
\end{align*}

For positive integers $k, m$ such that $2k \leq p$ and $1 \leq m \leq k$, we
define the $p \times p$ matrix $B(m, k) = (b_{ij})_{p\times{p}}$ such that:
$$
b_{ij} = \mathbf{1} \{i = m \text{ and } m+1 \leq j\leq 2k, \text{ or } 
            j = m \text{ and } m+1 \leq i \leq 2k\}
$$

Setting $k = n^\frac{1}{2\alpha + 1}$ and $a = k^{-\alpha-1)}$, we then define
a collection of $2^k$ covariance matrices:
$$
\mathcal{F}_{11} = \left\{
    \Sigma(\theta):
    \Sigma(\theta) = \mathbf{I}_{p\times p}
        + \tau a\sum_{m=1}^k \theta_m B(m, k),
    \quad\theta = (\theta_m) \in \{0, 1\}^k
\right\}
$$
where $0 < \tau < 2^{-\alpha - 1}M$.  We can interpret this parameter 
space intuitively.  First, we generate a set of $k$ matrices $B(m, k)$,
$m \in [k]$, in which $B(m, k)$ is nonzero only on certain elements
along row or column $m$.  Then, for each $\theta = \{0, 1\}^k$, we add
a number of $\tau a$ perturbations to the identity matrix; where
we add $\tau a B(m, k)$ if $\theta_m \neq 0$.  This gives us a set
of $\Sigma(\theta)$ of size $2^k$.

We now verify that $\mathcal{F}_{11} \subset \mathcal{F}(M, M_0)$.  Recall
the definition:
$$
\mathcal{F}(M_0, M) = \left\{
\Sigma: \max_j \sum_j \{|\sigma_{ij}|: |i-j| > k\} \leq M k^{-\alpha}
\text{ for all } k, \text{and} \lambda_{\text{max}}(\Sigma)\leq M_0
\right\}
$$
Let us consider the $\Sigma(\theta)$ where $\theta$ is a vector of 
$k$ ones and $\tau = 2^{-\alpha - 1}M$.  Then
$$
\Sigma(\theta) = \mathbf{I}_{p\times p}
    + M (k + 2)^{-\alpha -1}\sum_{m=1}^k B(m, k)
$$
We observe that $\sum_{m=1}^k B(m, k)$ is a matrix of ones on the
off-diagonals up to row and column $2k \leq p$.  Therefore, the worst-case
sum of entries more than $k$ away from the diagonal is:
\begin{align*}
    \frac{k}{2} M k^{-\alpha - 1} 2^{- \alpha - 1}
    &=  M k^{-\alpha} 2^{- \alpha - 2}  \\
    &\leq M k^{-\alpha}
\end{align*}
Without loss of generality the authors assume $M_0 > 1$ and $\rho > 1$; 
otherwise we may replace $\mathbf{I}_{p\times p}$ with $\varepsilon
\mathbf{I}_{ p \times p}$ with $\varepsilon \in (0, \min\{M_0, \rho\})$.  

[//]: # (TODO: Intepret this intuitively.)

To be continued!

