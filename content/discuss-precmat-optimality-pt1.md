Title: Discussion: Asymptotic Normality and Optimalities in Estimation of Large Gaussian Graphical Models
Date: 2017-02-09
Modified: 2017-02-09
Category: Statistics
Tags: math, stats, covariance, minimax, risk, precision, asymptotics, normality
Slug: discuss-precmat-optimality-pt1
Summary: Discussion of paper on asymptotic properties of an entrywise estimator for precision matrices.
<div style="display:none">
    $$
    \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    \newcommand{\RR}{\mathbf{R}}
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

The paper describes a technique for producing asymptotically efficient
entrywise estimators for precision matrices, under the Gaussian assumption.
They are able to accomplish this with a sparseness condition relative
to the sample size.  Intuitively, matrices that are less sparse require
more samples to achieve the parametric rate.  The parameter space also
makes an assumption the spectrum of the matrix (range of its singular
values).

In this initial post, some motivation for the problem is given, and then
we walk through some intuition, and then fully dissect the technique
used by the author, so that later posts can be fully dedicated to the
inference results and proofs of theorems.

## Motivation
The precision matrix is of particular interest in the multivariate Gaussian
setting due to the fact that the precision matrix admits, in the form of
its nonzero entries, an adjacency matrix for the graphical model associated
with conditional independence statements about the set of Gaussian random
variables involved.  This will be discussed in an upcoming blog post.

## Intuition
The authors produce entrywise elements of the precision matrix by exploiting
some neat properties about Gaussian vectors.

Consider a random matrix $X \in \RR^{n \times p}$, $n$ observations of $p$
variables, drawn from a $\Nn(0, \Sigma_{p\times p})$ distribution.  Rather
than estimating $\hat\Sigma$ directly and inverting it, we may make the
following observations.

Let's take a look at what happens if we take a pair of indices $A = \{i, j\}
\in [p]$, and regress the associated variables on all the others, $A^c = 
[p]\setminus\{i, j\}$.  We would have:
$$
x_A = \beta x_{A^c} + \epsilon
$$
where $\epsilon$ is a noise term, distributed normally with mean zero, and
which are independent of $A^c$.  Note that here

* $x_A \in \RR^{2}$
* $\beta \in \RR^{2\times (p-2)}$
* $x_{A^c} \in \RR^{p-2}$
* $\epsilon \in \RR^2$

We can multiply both sides by $x_{A^c}$, and then take the expectation on
both sides:
\begin{align*}
x_A &= \beta x_{A^c} + \epsilon \\
x_Ax_{A^c}^\top &= \beta x_{A^c}x_{A^c}^\top + \epsilon x_{A^c}^\top \\
\EE x_Ax_{A^c}^\top &= \beta \EE x_{A^c}x_{A^c}^\top    \\
\Sigma_{A, A^c} &=  \beta \Sigma_{A^c, A^c} \\
\Rightarrow
\beta &= \Sigma_{A, A^c}\Sigma_{A^c, A^c}^{-1}
\end{align*}

This immediately implies that given we observe $x_{A^c}$, the mean of $x_A$ is
$\Sigma_{A, A^c}\Sigma_{A^c, A^c}^{-1}x_{A^c}$.  But what about the variance?
Let's see how the noise is distributed:
\begin{align*}
\epsilon &= \Sigma_{A, A^c}\Sigma_{A^c, A^c}^{-1} x_{A^c} - x_A\\
\mathrm{var}(\epsilon) &= \Sigma_{A, A^c}\Sigma_{A^c, A^c}^{-1}\Sigma_{A^c,
    A^c}\Sigma_{A^c, A^c}^{-1}\Sigma_{A^c, A} + \Sigma_{A, A}\\ 
&=\Sigma_{A, A^c}\Sigma_{A^c, A^c}^{-1}\Sigma_{A^c, A} + \Sigma_{A, A}\\ 
\end{align*}

Then, assuming $x_{A^c}$ is observed, the variance of $x_A$ depends only
on the noise:
\begin{align*}
\mathrm{var}(x_A|x_{A^c})   &=  \mathrm{var}(\epsilon)  \\
&=\Sigma_{A, A^c}\Sigma_{A^c, A^c}^{-1}\Sigma_{A^c, A} + \Sigma_{A, A}\\ 
\end{align*}

We observe that the variance of $\epsilon$ is essentially the variance
of $x_A$ given $x_{A^c}$, which we denote $\Sigma_{A|A^c} = \Theta_{A, A}$.  By
using the fact that:
$$
    \Theta_{A, A} = \Omega_{A, A}^{-1}
$$
we can intuitively get estimates for blocks of the precision matrix at a time
by inverting the estimates for the conditional covariance matrix.

### Rewriting $\Sigma$ with $\Omega$
By crunching through some identities, we may rewrite:
$$
\beta = \Sigma_{A, A^c}\Sigma_{A^c, A^c}^{-1} 
    = -\Omega_{A, A}^{-1}\Omega_{A, A^c}
$$
which is similar to what the authors will use later on.  I can review the
salient identities in a future post.

## Methodology
Now that we've reviewed the intuition behind why this should work, let's dive
into the methodology of the paper.  Note that to be consistent with the
convention that observations are listed in a data matrix $X$ row-wise, most
of the vectors in the paper's exposition are row vectors; therefore, the
next few formulas will be transposes of those given in the previous section.

By analogy of the intuition previously given (the only differences are that
$\XX$ is now a matrix, and $\beta$ and $\epsilon$ are transposed), we have:
$$
\XX_A = \XX_{A^c}\beta + \epsilon_A
$$
where, adhering to the authors' usage of the precision matrix, and transposing
everything, we have:
$$
\beta = - \Omega_{A^c, A}\Omega_{A, A}
$$

Now, suppose we were only interested in estimating $\epsilon$ (as looking at
the noise is the key to estimating the entries of the precision matrix) and
we _knew_ $\beta$.  Then what would the maximum likelihood estimator of
$\Theta_{A, A}$ look like?  The authors denote oracle MLE of $\Theta_{A, A}$ as
$$
\Theta^{ora}_{A, A} = (\theta^{ora}_{kl})_{k, l \in A}
    = \frac{\epsilon_A^\top\epsilon_A}{n}
$$
and so the corresponding estimates of the precision matrix would be:
$$
\Omega^{ora}_{A, A} = (\omega^{ora}_{kl})_{k, l \in A}
    = \left(\Theta^{ora}_{A, A}\right)^{-1}
$$

But because $\beta$ is unknown, we must estimate $\beta$ and use its estimator
to then estimate $\epsilon_A$.  The authors introduce a scaled lasso regression
problem to do this.  For each $m \in A = \{i, j\}$, they perform the
optimization: 
$$
\left\{\hat\beta_m, \hat\theta^{1/2}_{mm}\right\}
=
\arg\min_{b\in\RR^{p-2},\\\sigma \in \RR^+}
\left\{
\frac{\norm{\XX_m - \XX_{A^c}b}^2}{2n\sigma}
+ \frac{\sigma}{2} 
+ \lambda\sum_{k\in A^c}\frac{\norm{\XX_k}}{\sqrt{n}}|b_k|
\right\}
$$

With these $\hat\beta$ estimates, we are able to define the residuals of the
scaled lasso regression:
$$
\hat\epsilon_A = \XX_A - \XX_{A^c}\hat\beta
$$
and consequently the estimate of the conditional covariance matrix:
$$
\hat\Theta_{A, A} = \frac{\hat\epsilon_A^\top\hat\epsilon_A}{n}
$$
and invert the entries of this estimator to get the estimates for $\Omega_{A,
A}$: 
$$
\hat\Omega_{A, A} = \hat\Theta_{A, A}^{-1}
$$

In my upcoming posts about this paper, I will review the inference results and
step through the associated proofs.
