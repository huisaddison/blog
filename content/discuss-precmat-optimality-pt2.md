Title: Discussion: Asymptotic Normality and Optimalities in Estimation of Large Gaussian Graphical Models, Part 2
Date: 2017-02-19
Modified: 2017-02-23
Category: Statistics
Tags: math, stats, covariance, minimax, risk, precision, asymptotics, normality
Slug: discuss-precmat-optimality-pt2
Summary: Discussion of paper on asymptotic properties of an entrywise estimator for precision matrices, continued.
<div style="display:none">
    $$
    \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    \newcommand{\RR}{\mathbf{R}}
    \newcommand{\DD}{\mathbf{D}}
    \newcommand{\PP}{\mathbf{P}}
    \newcommand{\EE}{\mathbf{E}}
    \newcommand{\XX}{\mathbf{X}}
    \newcommand{\Nn}{\mathcal{N}}
    \newcommand{\Nn}{\mathcal{N}}
    \DeclareMathOperator{\var}{var}
    \DeclareMathOperator{\diag}{diag}
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
$\ell_0$ norm; this special case is the parameter space used in **Theorem 1**.

Appealing to the graphical model representation of the multivariate normal
distribution, in which a nonzero entry $\omega_{ij}$ in the precision matrix
implies the existence of an edge between nodes $i$ and $j$, we may observe that
when $|\omega_{ij}|$ is zero or larger than $\lambda$, $s_\lambda$ is
equivalent to the maximum node degree of the graph (the degree of each node is
equivalent to the number of nonzero entries on each row, or column).  Finally,
we note that the spectrum (eigenvalues) of the matrix are bounded, a fact upon
which the later analysis relies.

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

First, we state a few conditions, which will be useful in our analysis of
**Theorem 2**.  When these conditions hold for certain fixed constant $C_0$,
$\varepsilon_\Omega \rightarrow 0$, and all $\delta \geq 1$, the asymptotic
normality and efficiency properties will hold, as we will see in the analysis
of Theorem 2.

The first condition is:
\begin{align}
\max_{A: A = \{i, j\}}
\PP\left\{
\norm{\XX_{A^c}\left(\hat\beta_{A^c, A} - \beta_{A^c, A}\right)}^2
\geq C_0 s \delta\log p
\right\}
\leq p^{-\delta + 1}\varepsilon_\Omega
\end{align}
Observing that $\XX_{A^c}\left(\hat\beta_{A^c, A} - \beta_{A^c, A}\right)$
is equivalent to $\norm{\epsilon_A - \hat\epsilon_A}^2$, we may interpret
this as a concentration bound on the deviation of the residual estimates
from the oracle residuals.

The next condition is:
\begin{align}
\max_{A:A = \{i, j\}}
\PP\left\{
\norm{
    \bar\DD^\frac{1}{2}_{A^c}
    \left(\hat\beta_{A^c, A} - \beta_{A^c, A}\right)
}_1 \geq C_0 s \sqrt{\delta\frac{\log p}{n}}
\right\}
\leq p^{-\delta+1}\varepsilon_\Omega
\end{align}
with $\bar \DD = \diag\left(\frac{\XX^\top\XX}{n}\right)$.  These two
statements are  essentially a risk bounds on the lasso estimator, which will be
discussed and proved in a future post.

The final condition is, for $\theta_{ii}^{ora} = \frac{\norm{\XX_i
- \XX{A^c}\beta_{A^c, i}}^2}{n}$,
\begin{align}
\max_{A: A = \{i, j\}}
\PP\left\{
\left|
\frac{\hat\theta_{ii}}{\theta_{ii}^{ora}} - 1
\right| \geq
C_0 s \delta \frac{\log p}{n}
\right\}
\leq p^{-\delta + 1}\varepsilon_\Omega
\end{align}
with a certain complexity measure $s$ of the precision matrix $\Omega$,
assuming the spectrum of $\Omega$ is bounded, and $n \geq \frac{(s\log
p)^2}{c_0}$ for a sufficiently small $c_0 > 0$.  

## Statement
**Theorem 2.** _Let $\hat\Theta_{A, A}$ and $\hat\Omega_{A, A}$ be estimators
defined in (2) and (3) respectively.  Let $\delta \geq 1$.  Suppose $s \leq
\frac{c_0n}{\log p}$ for a sufficiently small constant $c_0 > 0$._

1.  _Suppose that conditions (7), (8), (9) hold with $C_0$ and
    $\varepsilon_\Omega$.  Then
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
    estimators and $C_1$ is a positive constant depending only on
    $\{C_0, \max_{m\in A = \{i, j\}}\theta_{mm}\}$._
2.  _Let $\lambda = (1 + \varepsilon)\sqrt{\frac{2\delta\log p}{n}}$ with
    $\varepsilon > 0$ be the $\lambda$ parameter in the scaled lasso
    estimation problem, and let $\hat\beta_{A^c, A}$ be the scaled lasso 
    estimator, or the LSE after the scaled lasso selection.  Then (4),
    (5), and (6), and thus (7) and (8) hold for all $\Omega \in
    \mathcal{G}^*(M, s, \lambda)$ with a certain constant $C_0$ depending on
    $\{\varepsilon, c_0, M\} only and_
    \begin{align}
    \max_{\Omega \in \mathcal{G}^*(M, s, \lambda)} \varepsilon_\Omega = o(1)
    \end{align}
3.  _There exist constants $D_1$ and $\vartheta \in (0, \infty)$, and four
    marginally standard normal random variables $Z1, Z_{kl}$, where $kl = ii,
    ij, jj$, such that whenever $|Z_{kl}| \leq \vartheta\sqrt{n}$ for all $kl$,
    we have
    \begin{align}
    \left|\kappa_{ij}^{ora} - Z'\right|
    \leq
    \frac{D_1}{\sqrt{n}}\left(1 + Z_{ii}^2 + Z_{ij}^2 + Z_{jj}^2\right)
    \end{align}
    where $Z'$, which can be defined as a linear combination of $Z_{kl}$._

Intuitively, statement (1) says that there is a very low probability that
the maximum entrywise deviation of the actual estimator from the oracle MLE
is larger than a constant that we can control.  Statement (2) shows that
the conditions are met such that statement (1) holds.  Statement (3) says that
the rescaled oracle MLE behaves more or less asymptotically normally.

Once we show these statements about the estimator relative to oracle MLEs,
we will prove statements in *Theorem 3* relating the oracle MLEs to the true
parameter values, and by the triangle inequality, we will have bounds on
the distances between our estimates on the truth.

## Proof for Theorem 2(i)
The values of $\theta_{ii}, \theta_{jj}$ are uniformly bounded, which implies
that the desired concentration bound  (7) follows from (4) for
$\theta^{ora}_{ii}$ and $\theta^{ora}_{jj}$.  

Therefore, we only need to be concerned about bounding $\theta^{ora}_{ij}$.
Recall that we define $\bar \DD = \diag\left(\frac{\XX^\top \XX}{n}\right)$ and
that $\XX_{A^c}$ is independent of $\epsilon_A$.  First, we show the following.

**Claim.**

$$\left(\XX \bar \DD^{-\frac{1}{2}}\right)^\top_k
\frac{\epsilon_m}{n}\sim \Nn\left(0, \frac{\theta_{mm}}{n}\right)$$

_for all $m \in A$._

_Proof._  First, we observe that $\XX\bar\DD^{-\frac{1}{2}}$ is essentially
$\XX$ with its columns scaled to unit length in Euclidean norm.  The fact
that the mean of the distribution is zero follows from the fact that the
columns of $\XX$ are assumed to be centered.  To show the variance, we 
observe that we may express

\begin{align*}
\var\left(\left(\XX \bar \DD^{-\frac{1}{2}}\right)^\top_k\epsilon_m\right)
&=  \var\left(\sum_{i=1}^p \left(\XX\bar\DD^{-\frac{1}{2}}\right)_{ik}
    \epsilon_{im}\right)    \\
&=  \sum_{i=1}^p \left(\XX\bar\DD^{-\frac{1}{2}}\right)^2_{ik}
    \var\left(\epsilon_{im}\right)    \\
&=  \sum_{i=1}^p \left(\XX\bar\DD^{-\frac{1}{2}}\right)^2_{ik}
    \var\left(\epsilon_{1m}\right)    &&\text{(Symmetry.)}  \\
&=  \var\left(\epsilon_{1m}\right)  \\
&=  \EE \left[\epsilon_{1m}^2\right] &&\text{(Errors centered at zero.)}    \\
&=  n\theta_{mm}
\end{align*}

Dividing $\XX\bar\DD^{-1/2}$ by $n$ gives the desired
variance.  <div align="right"> &#8718; </div>

It then follows from the union bound and [Mill's Inequality](
{filename}mills-inequality.md) that:
$$
\PP\left\{
\norm{
\left(\XX\bar\DD^{-1/2}\right)^\top_{A^c} \frac{\epsilon_m}{n}
}_\infty > \sqrt{
2\delta\theta_{mm}n^{-1}\log p
}
\right\}
\leq \frac{p^{-\delta}(p-2)}{\sqrt{2\delta\log p}}
$$

Now, let's compare our covariance estimates to the oracle MLE:
\begin{align*}
\left| \hat\theta_{ij} - \theta_{ij}^{ora}\right|
    &=  \left|
        \frac{\hat\epsilon_i^\top\hat\epsilon_j}{n}
        -
        \frac{\epsilon_i^\top\epsilon_j}{n}
        \right|
\end{align*}
Recalling that
\begin{align*}
\hat\epsilon_A
    &= \XX_A - \XX_{A^c}\hat\beta_{A^c, A}  \\
\epsilon_A
    &= \XX_A - \XX_{A^c}\beta_{A^c, A}      \\
\Rightarrow
\hat\epsilon_A - \epsilon_A
    &=  \XX_{A^c}\left(\beta_{A^c, A} - \hat\beta_{A^c, A}\right)
\end{align*}
we have:
\begin{align*}
\left| \hat\theta_{ij} - \theta_{ij}^{ora}\right|
    &=  \frac{1}{n}\left|
        \hat\epsilon_i^\top\hat\epsilon_j
        -
        \epsilon_i^\top\epsilon_j
        \right| \\
    &=  \frac{1}{n}\left|
        \left(\epsilon_i 
        + (\hat\epsilon_i - \epsilon_i)
        \right)^\top
        \left(\epsilon_j
        + (\hat\epsilon_j - \epsilon_j)
        \right)
        -
        \epsilon_i^\top\epsilon_j
        \right| \\
    &=  \frac{1}{n}\left|
        \left(\epsilon_i 
        +
        \XX_{A^c}(\beta_i - \hat\beta_i)
        \right)^\top
        \left(\epsilon_j
        +
        \XX_{A^c}(\beta_j - \hat\beta_j)
        \right)
        -
        \epsilon_i^\top\epsilon_j
        \right| \\
    &=  \frac{1}{n}
        \Bigg|
        \epsilon_i^\top \epsilon_j
        + \left(\beta_i - \hat\beta_i\right)^\top \XX_{A^c}^\top \epsilon_j
        + \epsilon_i^\top \left(\beta_j - \hat\beta_j\right)\XX_{A^c}   \\
    &\qquad
        + \left(\beta_i - \hat\beta_i\right)^\top \XX_{A^c}^\top
            \XX_{A^c}\left(\beta_j - \hat\beta_j\right)
        -\epsilon_i^\top \epsilon_j
        \Bigg| \\
    &\leq  \frac{1}{n}\Bigg[
        \left|
        \left(\beta_i - \hat\beta_i\right)^\top \XX_{A^c}^\top \epsilon_j
        \right|
        + \left|
        \epsilon_i^\top \left(\beta_j - \hat\beta_j\right)\XX_{A^c}
        \right| \\
    &\qquad
        + \left|\left(\beta_i - \hat\beta_i\right)^\top \XX_{A^c}^\top
            \XX_{A^c}\left(\beta_j - \hat\beta_j\right)
        \right|\Bigg]  \\
    &=  \frac{1}{n}\Bigg[
        \left|
        \left(\beta_i - \hat\beta_i\right)^\top
        \bar\DD^{-1/2}\bar\DD^{1/2} \XX_{A^c}^\top \epsilon_j
        \right|
        + \left|
        \epsilon_i^\top \left(\beta_j - \hat\beta_j\right)
        \bar\DD^{-1/2}\bar\DD^{1/2}\XX_{A^c}
        \right| \\
    &\qquad
        + \left|\left(\beta_i - \hat\beta_i\right)^\top \XX_{A^c}^\top
            \XX_{A^c}\left(\beta_j - \hat\beta_j\right)
        \right|\Bigg]\\
    &\leq \frac{1}{n}\Bigg[
        \norm{
            \left(\XX\bar\DD^{-1/2}\right)_{A^c}^\top\epsilon_i
        }_\infty\norm{
            \bar\DD^{1/2}\left(\beta_j - \hat\beta_j\right)
        }_1
        +
        \norm{
            \left(\XX\bar\DD^{-1/2}\right)_{A^c}^\top\epsilon_j
        }_\infty\norm{
            \bar\DD^{1/2}\left(\beta_i - \hat\beta_i\right)
        }_1 \\
    &\qquad
        +
        \norm{
        \XX_{A^c}\left(\beta_i - \hat\beta_i\right)
        }\cdot\norm{
        \XX_{A^c}\left(\beta_j - \hat\beta_j\right)
        }
    \Bigg] \\
    &\leq 2\sqrt{2\delta\theta_{mm}n^{-1}\log p}C_0
            s\sqrt{\delta \frac{\log p}{n}}
            + \frac{C_0s\delta\log p}{n}\\
    &=  C_1 s\frac{\delta \log p}{n}
\end{align*}
with probability at least $1 - 2p^{-\delta + 1}\epsilon_\Omega
- 2p^{-\delta + 1}(2\log p)^{-1/2}$ by the union bound, implying (7).

Given that the spectrum of $\Theta_{A, A}$ is bounded, the functional
$\zeta_{kl}(\Theta_{A, A}) = \left(\Theta_{A, A}^{-1}\right)_{kl}$ is Lipschitz
in a neighborhood of $\Theta_{A, A}$ for $k, l \in A$, and thus the bound
on distances between the precision matrix estimates and the oracle MLE for
the precision matrix in (8) follows from (7).

## Proof for Theorem 2(ii)
The proof for part (ii), though fairly straightforward, depends on Theorem
10(i), Theorem 11(ii), and Proposition 1, and so we will return to this later.

This part of the Theorem essentially gives and proves conditions under which
conditions (4), (5), and (6) hold, which in turn imply (7) and (8) for all
$\Omega \in \Gg^*(M, s, \lambda)$ up to a constant.  This part of the Theorem
also establishes that $\varepsilon_\Omega$ is $o(1)$ for all $\Omega$ in the
parameter space, implying that it has a negligible impact on the concentration
bounds (7) and (8).

## Proof for Theorem 2(iii)
To prove the coupling inequality in (10), we first define a random vector
$\eta^{ora} = \left(\eta_{ii}^{ora}, \eta_{ij}^{ora}, \eta_{jj}^{ora}\right)$,
where
$$
\eta_{kl}^{ora} = \sqrt{n}\frac{\theta_{kl}^{ora} - \theta_{kl}}
{\theta_{kk}\theta_{ll} + \theta_{kl}^2}
$$
By the KMT inequality, for which the authors cite Mason and Zhou (2012) for the
one-dimensional case and Einmahl (1989) for the multidimensional case, there 
exist constants $D_0, \vartheta \in (0, \infty)$ and a random Gaussian vector
$Z = (Z_{ii}, Z_{ij}, Z_{jj}) \sim \Nn(0, \breve\Sigma)$ where $\breve\Sigma
= \cov(\eta^{ora})$, such that $|Z_{kl}| \leq \vartheta\sqrt{n}$ for all $kl$
implies
$$
\norm{\eta^{ora} - Z}_\infty \leq
\frac{D_0}{\sqrt{n}} \left(
1 + Z_{ii}^2 + Z_{ij}^2 + Z_{jj}^2
\right)
$$

Let us now define $\Theta = (\theta_{ii}, \theta_{ij}, \theta_{jj})$, consider
the function
$$
\omega_{ij}(\Theta) = -\frac{\theta_{ij}}{\theta_{ii}\theta_{jj} - \theta_{ij}^2}
$$
and take its Taylor expansion, which gives us:
\begin{align*}
    \omega_{ij}^{ora} - \omega_{ij}
    &=  \la \nabla\omega_{ij}(\Theta), \Theta^{ora} - \Theta\ra
        + \sum_{|\beta| = 2}R_\beta(\Theta^{ora})(\Theta-\Theta^{ora})^\beta
\end{align*}
where
\begin{align*}
    |\beta|     &\triangleq \sum_k \beta_k   \\
    x^\beta     &\triangleq \prod_k x_k^{\beta_k}    \\
    D^\beta f(x)&\triangleq \frac{\partial^{|\beta|} f}
        {
            \partial x_1^{\beta_1}
            \partial x_2^{\beta_2}
            \partial x_3^{\beta_3}
        }
\end{align*}
