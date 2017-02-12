Title: Discussion: Asymptotic Normality and Optimalities in Estimation of Large Gaussian Graphical Models, Part 2
Date: 2017-02-11
Modified: 2017-02-11
Category: Statistics
Tags: math, stats, covariance, minimax, risk, precision, asymptotics, normality
Slug: discuss-precmat-optimality-pt2
Summary: Discussion of paper on asymptotic properties of an entrywise estimator for precision matrices, continued.
<div style="display:none">
    $$
    \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    \newcommand{\RR}{\mathbf{R}}
    \newcommand{\PP}{\mathbf{P}}
    \newcommand{\EE}{\mathbf{E}}
    \newcommand{\XX}{\mathbf{X}}
    \newcommand{\Nn}{\mathcal{N}}
    \{\Nn}{\mathcal{N}}
    \DeclareMathOperator{\var}{var}
    $$
</div>

## Introduction
As part of project I've been working on, I'm reading [Asymptotic Normality and
Optimalities in Estimation of Large Gaussian Graphical Models](
https://arxiv.org/abs/1309.6024), a paper by Ren, Sun, Zhang, and Zhou.

Motivation, problem setup, and other preliminaries are addressed in a previous
[blog post]({filename}discuss-precmat-optimality-pt1.md).  In this post, we
will begin discussing the statistical inference results of the paper by
walking through the proof for **Theorem 2**, which places bounds on the
distribution of the estimates and makes the results in **Theorem 3** possible.
The results of **Theorem 1** are simply a special case of **Theorem 2**.

## Preliminaries
First, we recall the definitions of the estimators for the sub-blocks of the
precision and conditional covariance matrices.  We begin with the scaled
lasso regression parameter estimates:  
\begin{align}
\left\{\hat\beta_m, \hat\theta^{1/2}_{mm}\right\}
=
\arg\min_{b\in\RR^{p-2},\\\sigma \in \RR^+}
\left\{
\frac{\norm{\XX_m - \XX_{A^c}b}^2}{2n\sigma}
+ \frac{\sigma}{2} 
+ \lambda\sum_{k\in A^c}\frac{\norm{\XX_k}}{\sqrt{n}}|b_k|
\right\}
\end{align}
which yield $\hat\epsilon_A$ as the residual estimates from regression $\XX_A$
on $\XX_{A^c}$.  Then, we define:
\begin{align}
\hat\Theta_{A, A}   &=  \frac{\hat\epsilon_A^\top\epsilon_A}{n} \\
\hat\Omega_{A, A}   &=  \hat\Theta_{A, A}^{-1}
\end{align}

We now define the parameter space considered.  For $\lambda > 0$, we
defined capped-$\ell_1$ balls as follows:
$$
\mathcal{G}^* = \left\{
\Omega: s_\lambda(\Omega) \leq s, M^{-1}
    \leq \lambda_\min(\Omega)
    \leq \lambda_\max(\Omega)
    \leq M
\right\}
$$
where
$$
s_\lambda = s_\lambda(\Omega) = \max_j\sum_{i\neq j}
\min\left\{1, \frac{|\omega_{ij}|}{\lambda}\right\}
$$
for $\Omega = (\omega_{ij})_1\leq i, j\leq p$.  The authors take $\lambda$
on the order $\sqrt\frac{\log p}{n}$ in this paper.

Intuitively, this parameter space is a mix of the $\ell_0$ norm, which
measures sparsity, and the $\ell_1$ norm, which is the sum of absolute
values, imposed on each row.  For $\lambda$ very small, we recover the pure
$\ell_0$ norm; this special case is **Theorem 1**.

Appealing the graphical model representation of the multivariate normal
distribution with nonzero entries of the precision matrix giving the edges,
we may observe that when $|\omega_{ij}|$ is zero or larger than $\lambda$, 
$s_\lambda$ is equivalent to the maximum node degree of the graph (the
degree of each node is equivalent to the number of nonzero entries on
each row, or column).  Finally, we note that the spectrum (eigenvalues) of
the matrix are bounded, a fact upon which the later analysis relies.

The authors then prove a theorem that gives an error bound on the estimates.
Their approach is to:

1.  Compare the estimates to the oracle MLE, giving a concentration bound
    on the distances between them.
2.  Show that
    $$
    \kappa_{ij}^{ora} = \sqrt{n}
    \frac{\omega_{ij}^{ora} - \omega_{ij}}
    {\sqrt{\omega_{ii}\omega_{jj} + \omega_{ij}^2}}
    $$
    is asymptotically standard normal, which implies the oracle MLE is
    asymptotically normal with mean $\omega_{ij}$ and variance
    $n^{-1}\sqrt{\omega_{ii}\omega_{jj} + \omega_{ij}^2}$.

By coupling the actual estimator to the oracle MLE and then proving nice
properties for the oracle MLE, we can work towards nice properties for
the actual estimator.

## Statement
**Theorem 2.** _Let $\hat\Theta_{A, A}$ and $\hat\Omega_{A, A}$ be estimators
defined in (2) and (3) respectively, and $\lambda = (1 + \varepsilon)
\sqrt\frac{2\delta\log p}{n}$ for any $\delta \geq 1$ and $\varepsilon > 0$
in equation (1)._

1.  _Suppose $s \leq \frac{c_0n}{\log p}$ for a sufficiently small constant
    $c_0 > 0$.  We have
    \begin{align}
    \max_{G^*(M, s, \lambda)}\max_{A:A=\{i, j\}}
    \PP\left\{
        \norm{\hat\Theta_{A, A} - \Theta_{A, A}^{ora}}_\infty
        > C_1 s \frac{\log p}{n}
    \right\}\leq o\left(p^{-\delta + 1}\right)
    \end{align}
    and
    \begin{align}
    \max_{G^*(M, s, \lambda)}\max_{A:A=\{i, j\}}
    \PP\left\{
        \norm{\hat\Omega_{A, A} - \Omega_{A, A}^{ora}}_\infty
        > C_1' s \frac{\log p}{n}
    \right\}\leq o\left(p^{-\delta + 1}\right)
    \end{align}
    where $\Theta^{ora}_{A, A}$ and $\Omega^{ora}_{A, A}$ are the oracle
    estimators and $C_1$ and $C_1'$ are positive constants depending on
    $\{\varepsilon, c_0, M\}$ only._
2.  _There exist constants $D_1$ and $\vartheta \in (0, \infty)$, and three
    marginally standard normal random variables $Z_{kl}$, where $kl = ii, ij,
    jj$, such that whenever $|Z_{kl}| \leq \vartheta\sqrt{n}$ for all $kl$,
    we have
    \begin{align}
    \left|\kappa_{ij}^{ora} - Z'\right|
    \leq
    \frac{D_1}{\sqrt{n}}\left(1 + Z_{ii}^2 + Z_{ij}^2 + Z_{jj}^2\right)
    \end{align}
    where $Z' \sim \Nn(0, 1)$, which can be defined as a linear combination
    of $Z_{kl}$._

Intuitively, statement (1) says that there is a very low probability that
the maximum entrywise deviation of the actual estimator from the oracle MLE
is larger than a constant that we can control.  Statement (2) says that
the rescaled oracle MLE behaves more or less asymptotically normally.

[//]:   #(TODO: Talk to Harry to get more intuition about these statements, update.)

## Proof for Theorem 2.1
First, we present the following lemma:

**Lemma 2.** _Let $\lambda = (1 + \varepsilon)\sqrt\frac{2\delta\log p}{n}$ for
any $\delta \geq 1$ and $\varepsilon > 0$.  Define the event $E_m$ as follows:
\begin{align}
\left|
\hat\theta_{mm} - \theta^{ora}_{mm} \right|
    &\leq   C_1'\lambda^2s  \\
\norm{\beta_m - \hat\beta_m}_1
    &\leq C_2'\lambda s   \\
n^{-1}\norm{\XX_{A^c}\left(\beta_m - \hat\beta_m\right)}^2
    &\leq C_3'\lambda^2 s   \\
\norm{n^{-1}\XX_{A^c}^\top\epsilon_m}_\infty
    &\leq C_4'\lambda
\end{align}
for $m \in \{i, j\}$ and some constants $C_k', 1 \leq k \leq 4$.  Under
the assumptions of Theorem 2, we have:_
$$
\PP\left\{E_m^c\right\} \leq o\left(p^{-\delta + 1}\right)
$$

The proof of **Lemma 2** is deferred to a future post.

[//]:   #(TODO: Read and understand the proof for Lemma 2.)


First, we consider $\theta_{ii}^{ora}$ and $\theta_{jj}^{ora}$.  By the
definition of $\lambda$ and Equation (7) of **Lemma 2**, the concentration
bound on the deviations in (4) holds for $\theta_{ii}^{ora},
\theta_{ij}^{ora}$.  We then consider the event $E_i \cap E_j$:

\begin{align*}
\left|\hat\theta_{ij} - \theta^{ora}_{ij}\right|
&=  \left|
    \frac{\hat\epsilon_i^\top\hat\epsilon_j}{n}
    - 
    \frac{\epsilon_i^\top\epsilon_j}{n}
    \right|
\end{align*}

(To be continued!)

## Proof for Theorem 2.2
