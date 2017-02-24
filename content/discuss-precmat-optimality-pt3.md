Title: Discussion: Asymptotic Normality and Optimalities in Estimation of Large Gaussian Graphical Models, Part 3
Date: 2017-02-23
Modified: 2017-02-23
Category: Statistics
Tags: math, stats, covariance, minimax, risk, precision, asymptotics, normality
Slug: discuss-precmat-optimality-pt3
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
    \newcommand{\Gg}{\mathcal{G}}
    \newcommand{\la}{\langle}
    \newcommand{\ra}{\rangle}
    \DeclareMathOperator{\var}{var}
    \DeclareMathOperator{\diag}{diag}
    \DeclareMathOperator{\cov}{cov}
    $$
</div>

## Introduction
The post reviews completes the proof for the risk upper bound on the Asymptotic
Normal Thresholding procedure introduced in a [paper by Zhou et al](
https://arxiv.org/abs/1309.6024).  The main object of interest is **Theorem
3**, which bounds the distance from the oracle MLE from the truth in
probability.

Previous posts give a [general overview of the paper's results](
{filename}discuss-precmat-optimality-pt1.md) and a [walkthrough of **Theorem
2**]({filename}discuss-precmat-optimality-pt2.md), which bounds the distance
between the proposed estimates and the oracle MLEs in probability.

## Statement
**Theorem 3.** _Let $\hat\Omega_{A, A}$ be the estimator of $\Omega_{A, A}$
defined in (2) in [Post 2]({filename}discuss-precmat-optimality-pt1.md) with
the components of $\hat\varepsilon_A$ being the estimated residuals of the
scaled lasso solution or the LSE of the selected model.  Set $\lambda = 
(1 + \varepsilon)\sqrt{\frac{2\delta\log p}{n}}$ in the scaled lasso 
problem with certain $\delta \geq 1$ and $\varepsilon > 0$.  Suppose $s \leq
c_0\frac{n}{\log p}$ for a sufficiently small constant $c_0 > 0$.  For any
small constant $\varepsilon_0 > 0$, there exists a constant $C_2 = 
C_2(\varepsilon_0, \varepsilon, c_0, M)$ such that
\begin{align}
\max_{\Omega\in\Gg^*(M, s, \lambda)}\max_{1 \leq i \leq j \leq p}
\PP\left\{
|\hat\omega_{ij} - \omega_{ij}| > C_2 \max\left\{
s\frac{\log p}{n}, \sqrt\frac{1}{n}
\right\}
\right\} \leq \varepsilon_0
\end{align}
Morever, there exists a constant $C_3 = C_3(\delta, \varepsilon, c_0, M)$ such
that
\begin{align}
\max_{\Omega\in\Gg^*(M, s, \lambda)}
\PP\left\{
\norm{\hat\Omega - \Omega}_\infty > C_3 \max\left\{
s\frac{\log p}{n}, \sqrt\frac{\log p}{n}
\right\}
\right\} \leq o(p^{-\delta + 3})
\end{align}
Furthermore, $\hat\omega_{ij}$ is asymptotically efficient with a consistent
variance estimate
\begin{align}
\sqrt{nF_{ij}}(\hat\omega_{ij} - \omega_{ij})
\stackrel{D}{\rightarrow}\Nn(0, 1),
\qquad \frac{\hat F_{ij}}{F_{ij}} \rightarrow 1
\end{align}
uniformly for all $i, j$ and $\Omega \in \Gg^*(M, s, \lambda)$, provided that
$s = o\left(\frac{\sqrt{n}}{\log p}\right)$, where_
\begin{align}
F_{ij} = (\omega_{ii}\omega_{jj}- \omega_{ij}^2)^{-1},
\qquad
\hat F_{ij} = (\hat\omega_{ii}\hat\omega_{jj}- \hat\omega_{ij}^2)^{-1}
\end{align}
