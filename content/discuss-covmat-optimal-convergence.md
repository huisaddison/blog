Title: Discussion: Optimal Rates of Convergence for Covariance Matrix Estimation
Date: 2017-01-26
Category: Statistics
Tags: math, stats, covariance, minimax, risk
Slug: discuss-covmat-optimal-convergence
Summary: Discussion of paper on minimax estimation of covariance matrices.

For the past couple of days, I have been reading [Optimal Rates of Convergence
for Covariance Matrix Estimation](https://arxiv.org/abs/1010.3866).  In this
post, we will discuss the results from the paper and walk through steps of the
proofs.

The authors in this paper establish optimal rates of convergence for estimating
the covariance matrix under the operator norm (which, in Euclidean distance,
coincides with the spectral norm) and the Frobenius norm; i.e., they provide
lower bounds on the max risk and show that the minimax upper bound of their
tapering estimators achieves this rate.

First, they construct the following parameter space of $p\times{p}$ covariance
matrices:
$$
\mathcal{F}(M_0, M) = \left\{
\Sigma: \max_j \sum_j \{|\sigma_{ij}|: |i-j| > k\} \leq M k^{-\alpha}
\text{ for all } k, \text{and} \lambda_{\text{max}}(\Sigma)\leq M_0
\right\}
$$
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

To estimate $\Sigma$ by tapering the maximum likelihood estimator:
$$
\hat\Sigma = \hat\Sigma_k = \left(w_{ij}\sigma^*_{ij}\right)_{p\times p}
$$
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
each side at $k/2$.  

(To be continued!)
