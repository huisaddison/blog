Title: Discussion: Optimal Rates of Convergence for Covariance Matrix Estimation, Part 1
Date: 2017-01-26
Modified: 2017-02-01
Category: Statistics
Tags: math, stats, covariance, minimax, risk
Slug: discuss-covmat-optimal-convergence-pt1
Summary: Discussion of paper on minimax estimation of covariance matrices.

For the past couple of days, I have been reading [Optimal Rates of Convergence
for Covariance Matrix Estimation](https://arxiv.org/abs/1010.3866).  To aid my
digestion of this paper, I will be writing about it on this blog.  This post
covers the authors' justification of **Theorem 2**.  I highly recommend
referring to the paper; I will try to fill in gaps of the proof, but must leave
out some details for brevity.

## Background
The authors in this paper establish optimal rates of convergence for estimating
the covariance matrix under the operator norm (which, in Euclidean distance,
coincides with the spectral norm) and the Frobenius norm; i.e., they provide
lower bounds on the max risk and show that the minimax upper bound of their
tapering estimators achieves this rate.

First, they construct the following parameter space of $p\times{p}$ covariance
matrices:
\begin{equation}\label{eq:parameter-space}
\mathcal{F}(M_0, M) = \left\{
\Sigma: \max_j \sum_j \{|\sigma_{ij}|: |i-j| > k\} \leq M k^{-\alpha}
\text{ for all } k, \text{and} \lambda_{\text{max}}(\Sigma)\leq M_0
\right\}
\end{equation}
This definition states that for the covariance matrices in this space, the sum
of the absolute values of the entries $k$ away from the diagonal decay with
$k^{-\alpha}$, where $\alpha$ is a smoothing parameter.  The bound
will play an important proof in the later proofs.

**Theorem.** The minimax risk of estimating the covariance matrix $\Sigma$
over the class $\mathcal{P}_\alpha$ given above satisfies
$$
\inf_{\hat\Sigma}\sup_{\mathcal{P}_\alpha}
\mathbf{E}\Vert \hat\Sigma - \Sigma\Vert^2 \asymp
\min\left\{n^{-\frac{2\alpha}{2\alpha+1}}+\frac{\log p}{n}, \frac{p}{n}\right\}
$$
where $\mathcal{P}_\alpha = \mathcal{P}_\alpha(M_0, M, \rho)$ is the set of all
distributions of $X$ that satisfies both Equations \eqref{eq:parameter-space}
and \eqref{eq:subgaussian}

$\Sigma$ is then estimated by tapering the maximum likelihood estimator:
\begin{equation}\label{eq:tapering-estimator}
\hat\Sigma = \hat\Sigma_k = \left(w_{ij}\sigma^*_{ij}\right)_{p\times p}
\end{equation}
where $\sigma^*_{ij}$ are the entries in the maximum likelihood estimator
$\Sigma^*$ and the weights are given by:
\begin{align*}
w_{ij} = \begin{cases}
    1   &\text{ for } |i - j| \leq \frac{k}{2}   \\
    2 - 2\frac{|i-j|}{k} &\text{ for } \frac{k}{2} \leq |i - j| \leq k  \\
    0   &\text{ otherwise }
\end{cases}
\end{align*}

The bandwidth is $k$ on either side along the diagonal; shrinkage begins on
each side at $k/2$.  As a technical note, such a tapering estimator may be
rewritten as a sum of block matrices along the diagonal; this is used in the
proofs to attain concentration bounds via random matrix theory.

## Minimax Upper Bound under the Operator Norm
The authors assume that the $X_i$'s are subgaussian; that is, there exists
$\rho > 0$ such that:
\begin{equation}\label{eq:subgaussian}
\mathbf{P}\{|v^T(X_1 - \mathbf{E} X_1)| > t\} \leq \exp\left\{
-\frac{t^2\rho}{2}\right\}
\end{equation}
for all $t > 0$ and $\Vert{v}\Vert_2=1$.

**Theorem.** _The tapering estimator $\hat\Sigma_k$ defined in Equation
\eqref{eq:tapering-estimator} with $p \geq n^{\frac{1}{2\alpha+1}}$ satisfies_
\begin{equation}
\sup_{\mathcal{P}_\alpha}\mathbf{E}\Vert \hat\Sigma_k - \Sigma\Vert^2 \leq C
\frac{k+\log p}{n} + Ck^{-2\alpha}
\end{equation}
_for $k = O(n), \log p = O(n)$ and some constant $C > 0$._

To prove this, the authors assume $\mu = 0$ and analyze
$$
\tilde \Sigma = \frac{1}{n}\sum_{i=1}^n X_l X_l^\top
$$
rather than the maximum likelihood estimator
$$
\Sigma^* = \frac{1}{n}\sum_{i=1}^n X_lX_l^\top - \bar X \bar X^\top
$$
as $\bar X \bar X^\top$ is a higher order term and can be neglected in the
analysis of the rate.  As such, we defined a new tapering estimator:
$$
\breve\Sigma = \left(\breve\sigma_{ij}\right)_{1 \leq i,j \leq p}
= \left(w_{ij}\tilde\sigma_{ij}\right)_{1 \leq i,j \leq p}
$$

We observe that we may bound the risk from above by the bias and variance:
\begin{align*}
    \mathbf{E} \Vert \breve\Sigma - \Sigma\Vert^2
    &\leq  2\mathbf{E}\Vert \breve\Sigma - \mathbf{E} \breve\Sigma\Vert^2
        + 2\Vert\mathbf{E}\breve\Sigma- \Sigma\Vert^2
\end{align*}
This inequality can easily be derived by thinking about how to bound $(a +
b)^2$.

### Bias
To prove the bound on the bias, the authors note that the operator norm of a
symmetric matrix is upper bounded by its $\ell_1$ norm.  Therefore, we have:
\begin{equation}\label{eq:operator-upperbound-bias}
\Vert \mathbf{E}\breve\Sigma - \Sigma\Vert^2
\leq \left[\max_{i=1, \dotsc, p} \sum_{j: |i-j| > k} |\sigma_{ij}|\right]^2
\leq M^2 k^{-2\alpha}
\end{equation}
This inequality holds by construction (see Equation \eqref{eq:parameter-space}).

### Variance
The authors rely on random matrix theory to bound the variance.  Of particular
note is

**Lemma 2.** _There exists a submatrix $M_l^{(m)}$:_
$$
\Vert \breve\Sigma - \mathbf{E}\breve\Sigma\Vert \leq 3 N^{(m)}
= \Vert M_l^{(m)} - \mathbf{E} M_l^{(m)}\Vert
$$

They also provide a concentration bound on the operator norm of this submatrix.

**Lemma 3.** _There is a constant $\rho_1 > 0$ such that_
$$
\mathbf{P}\left\{N^{(m)} > x\right\} \leq 2p5^m\exp(-nx^2\rho_1)
$$
_for all $0 < x < \rho_1$ and $1-m \leq 1 \leq p$._

Therefore, by Lemma 2, we have:
\begin{align*}
\mathbf{E}\Vert \breve\Sigma - \mathbf{E}\breve\Sigma\Vert^2
&\leq 9 \mathbf{E}\left(N^{(m)}\right)^2                    \\
&= 9 \mathbf{E}\left(N^{(m)}\right)^2[I(N^{(m)} \leq x) + I(N^{(m)} > x)]   \\
&\leq 9 [x^2 +\mathbf{E} \left(N^{(m)}\right)^2I(N^{(m)} > x)]
\end{align*}

The next few steps are quite tricky to understand.

I believe that we can recall the definition of $N^{(m)}$:
$$
N^{(m)} = \Vert M_l^{(m)} - \mathbf{E} M_l^{(m)}\Vert
$$
for some $l$.  We note that we may bound the operator norm of a submatrix by
the operator norm of the full matrix by recalling the definition of an
operator norm
$$
\Vert Av \Vert_{\mathrm{op}} = \sup\{\Vert{Av}\Vert: v \in V, \Vert v \Vert = 1\}
$$
and seeing that for any vector $v$, $\Vert A_\mathrm{sub} v \Vert \leq \Vert
Av\Vert$.  Therefore, we may establish that:
$$
\Vert M_l^{(m)} - \mathbf{E} M_l^{(m)}\Vert
\leq
\Vert \breve\Sigma - \mathbf{E}\breve\Sigma\Vert
$$
at this point, we may split the norm by the triangle inequality:
$$
\Vert \breve\Sigma - \mathbf{E}\breve\Sigma\Vert
\leq
\Vert \breve\Sigma\Vert + \Vert\mathbf{E}\breve\Sigma\Vert
$$
and bound the operator norm by the Frobenius norm:
$$
\Vert \breve\Sigma\Vert + \Vert\mathbf{E}\breve\Sigma\Vert
\leq \Vert\breve\Sigma\Vert_F + C
$$
Finally, we may apply Cauchy-Schwarz:
\begin{align*}
\mathbf{E}\Vert \breve\Sigma - \mathbf{E}\breve\Sigma\Vert^2
&\leq C_1 \left[x^2 + \mathbf{E}\left(\Vert\breve\Sigma\Vert_F + C\right)
I(N^{(m)} > x)\right]   \\
&\leq C_1 \left[x^2 + \sqrt{\mathbf{E}\left(\Vert\breve\Sigma\Vert_F +
C\right)^4} \sqrt{\mathbf{P}(N^{(m)} > x)}\right]
\end{align*}

We now bound the $\sqrt{\mathbf{P}(N^{(m)} > x)}$ term.  By setting
$x=4\sqrt{\frac{\log{p}+m}{n\rho_1}}$, and recalling **Lemma 3**, we have:
\begin{align*}
\sqrt{\mathbf{P}(N^{(m)} > x)}
&\leq \sqrt{2p5^m \exp\{-nx^2\rho_1\}}  \\
&\leq \sqrt{2p5^m \exp\{-16\log p + -16m\}}  \\
&\leq \sqrt{2p5^m \cdot p^{-16}\exp\{-16m\}}  \\
\end{align*}

We are able to bound the Frobenius norm:
$$
\sqrt{\mathbf{E}\left(\Vert\breve\Sigma\Vert_F + C\right)^4} \leq p^2
$$
by observing that the squared Frobenius norm decouples the entries of the
matrix, and we are able to bound each of the $p^2$ entries by a constant,
which hides among the other constants.

With all these pieces together, we may conclude:
\begin{align*}
\mathbf{E}\Vert \breve\Sigma - \mathbf{E}\breve\Sigma\Vert^2
&\leq C\left[
    \frac{\log p + m}{n} + \underbrace{
        p^2\cdot(p5^m)^\frac{1}{2}\cdot p^{-8}
        \exp\{-8m\}
    }_\text{Lower Order Term}
    \right] \\
&\leq C_1 \left(\frac{\log p + m}{n}\right)
\end{align*}
where $C$ is an evolving constant.

### Putting It All Together
Having walked through the steps of bounding the bias:
$$
\Vert \mathbf{E}\breve\Sigma - \Sigma\Vert^2
\leq M^2 k^{-2\alpha} = Ck^{-2\alpha}
$$

and variance:
$$
\mathbf{E}\Vert \breve\Sigma - \mathbf{E}\breve\Sigma\Vert^2
\leq C \left(\frac{\log p + m}{n}\right)
$$
for arbitrary $\Sigma$, we can put these terms together to establish the bound
on the worse case risk of the tapering estimator (**Theorem 2**):
\begin{equation}
\sup_{\mathcal{P}_\alpha}\mathbf{E}\Vert \hat\Sigma_k - \Sigma\Vert^2 \leq C
\frac{k+\log p}{n} + Ck^{-2\alpha}
\end{equation}
for $k = O(n), \log p = O(n)$ and some constant $C > 0$.

