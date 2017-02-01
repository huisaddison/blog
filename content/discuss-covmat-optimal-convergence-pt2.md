Title: Discussion: Optimal Rates of Convergence for Covariance Matrix Estimation, Part 2
Date: 2017-01-30
Modified: 2017-01-31
Category: Statistics
Tags: math, stats, covariance, minimax, risk
Slug: discuss-covmat-optimal-convergence-pt2
Summary: Discussion of paper on minimax estimation of covariance matrices, continued.
$$
    \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
$$

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
\mathbf{I}_{ p \times p}$ with $\varepsilon \in (0, \min\{M_0, \rho\})$.  It
follows that $\mathcal{F_{11}} \subset \mathcal{F}(M, M_0)$.

We now construct a second parameter space $\mathcal{F}_{12} \subset
\mathcal{F}$ as follows:
$$
\mathcal{F}_{12} = \left\{
\Sigma_m: \Sigma_m \mathbf{I}_{p\times p} + \left(
\sqrt{\frac{\tau}{n}}\mathbf{1}\{i = j = m\}\right)_{p \times p},
0 \leq m \leq p_1
\right\}
$$
where $p_1 = \min\{p, e^\frac{n}{2}\}$ and $0 < \tau < \min\{(M_0 - 1)^2,
(\rho - 1)^2, 1\}$.  Because the $\Sigma_m$ is this parameter are diagonal
matrices, the bandability condition is satisfied trivially, and because the
greatest diagonal entry is $1 + \sqrt{\frac{\tau}{n}}$, the condition that
the spectral norm be less than $M_0$ is easily verified.

Denote $\mathcal{F}_1 = \mathcal{F}_{11} \cup \mathcal{F}_{12}$.  Observe
that $\mathcal{F}_1 \subset \mathcal{F}_\alpha(M_0, M)$

We will now prove that:
$$
\inf_{\hat\Sigma}\sup_{\mathcal{F}_{11}}
\mathbf{E} \norm{\hat\Sigma - \Sigma}^2 \geq cn^{-\frac{2\alpha}{2\alpha+1}}
$$
and
$$
\inf_{\hat\Sigma}\sup_{\mathcal{F}_{12}}
\mathbf{E} \norm{\hat\Sigma - \Sigma}^2 \geq c\frac{\log p}{n}
$$
for some constant $c > 0$.  Taken together, we have:
$$
\inf_{\hat\Sigma}\sup_{\mathcal{F}_1}
\mathbf{E} \norm{\hat\Sigma - \Sigma}^2
\geq \frac{c}{2}\left(n^{-\frac{2\alpha}{2\alpha+1}} + \frac{\log p}{n}\right)
$$
which proves **Theorem 3**.

## Lower Bound by Assouad's Lemma
Now that we have our machinery set up, we can move on the meat of our proof.
The goal of this section is to establish:

$$
\inf_{\hat\Sigma}\sup_{\mathcal{F}_{11}}
\mathbf{E} \norm{\hat\Sigma - \Sigma}^2 \geq cn^{-\frac{2\alpha}{2\alpha+1}}
$$

Suppose we are estimating an arbitrary quantity $\psi(\theta)$ within a
metric space with metric $d$, over a set of parameters $\Theta = \{0, 1\}^k$.
We denote the Hamming distance on $\{0, 1\}^k$ by $H(\theta, \theta') = \sum_{i
= 1}^k |\theta_i - \theta_i'|$.  We also define the total variation affinity
between two probability measures $P$ and $Q$ with densities $p, q$ with respect
to measure $\mu$ by $\norm{P \wedge Q} = \int \min\{p, q\} d\mu$.  Under these
assumptions, Assoud's Lemma gives the following lower bound on the maximum risk
of estimating $\psi(\theta)$:

**Lemma (Assouad).**  _Let $\Theta = \{0, 1\}^k$ and let $T$ be an estimator
based on an observation from a distribution in the collection $\{P_\theta,
\theta \in \Theta\}$.  Then for all $s > 0$:_
$$
\max_{\theta\in\Theta} 2^s \mathbf{E}_\theta d^s(T, \psi(\theta))
\geq \min_{H(\theta, \theta') \geq 1}
\frac{
    d^s(\psi(\theta), \psi(\theta'))
}{
    H(\theta, \theta')
}\cdot\frac{k}{2}\cdot
\min_{H(\theta, \theta')}\norm{\mathbf{P}_\theta \wedge\mathbf{P}_{\theta'}}
$$

The authors give a natural interpretation of Assouad's Lemma in terms of
multiple comparisons:

[//]: # (Why are there only k comparisons?)

1.  The first factor is the minimum cost of making a mistake per comparison;
    that is, it is a lower bound on the distance between the distance between
    two parameters in the parameter space.
2.  The last factor (total variation affinity) measures the overlap between
    the two probability measures indexed by $\theta$ and $\theta'$;
    intuitively, it gives a lower bound on the total probability of making
    type I and type II errors for each comparison.
3.  $\frac{k}{2}$ is the expected number of mistakes made when
    $\mathbf{P}_\theta$ and $\mathbf{P}_{\theta'}$ are indistinguishable
    from each other when $H(\theta, \theta') = 1$.

Assouad's Lemma can be unpacked as an extension of Le Cam's Method, and will
be discussed in an upcoming post, with guidance from [Yu's _Assoud, Fano,
and Le Cam_](https://www.stat.berkeley.edu/~binyu/ps/LeCam.pdf).

Suppose we draw $\mathbf{X}_1, \dotsc, \mathbf{X}_n \stackrel{\text{iid}}{\sim}
\mathcal{N}(0, \Sigma(\theta))$ with $\Sigma(\theta) \in \mathcal{F}_{11}$.  We
denote the joint distribution for these random vectors by $\mathbf{P}_\theta$.
The authors give two lemmas to complete the proof:

**Lemma 5.**  _Let $\Sigma(\theta)\in\mathcal{F}_{11}$.  Then for some constant
$c > 0$:_
$$
    \min_{H(\theta, \theta') \geq 1}
    \frac{
    \norm{\Sigma(\theta) - \Sigma(\theta')}^2
    }{
    H(\theta, \theta')
    } \geq cka^2
$$

**Lemma 6.**  _Suppose we draw $\mathbf{X}_1, \dotsc, \mathbf{X}_n
\stackrel{\text{iid}}{\sim} \mathcal{N}(0, \Sigma(\theta))$ with
$\Sigma(\theta) \in \mathcal{F}_{11}$.  Denote the joint distribution
by $\mathbf{P}_\theta$.  Then for some constant $c > 0$:_
$$
    \min_{H(\theta, \theta') = 1}
    \norm{\mathbf{P}_\theta \wedge \mathbf{P}_{\theta'}} > c
$$
By Lemmata 5 and 6, and taking $k = n^\frac{1}{2\alpha+1}$, we have the
desired bound:
$$
\sup_{\mathcal{F}_{11}}
\mathbf{E} \norm{\hat\Sigma - \Sigma}^2
\geq \frac{c^2}{2}k^2a^2
\geq c_1n^{-\frac{2\alpha}{2\alpha+1}}
$$
for some $c_1 > 0$.  As the bound is for an arbitrary $\hat\Sigma$, it
follows that:
$$
\inf_{\hat\Sigma}\sup_{\mathcal{F}_{11}}
\mathbf{E} \norm{\hat\Sigma - \Sigma}^2 \geq cn^{-\frac{2\alpha}{2\alpha+1}}
$$

## Lower Bound by Le Cam's Method

To be continued!

