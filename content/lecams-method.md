Title: Le Cam's Method
Date: 2017-02-05
Modified: 2017-02-05
Category: Statistics
Tags: math, stats, minimax, risk, le-cam
Slug: lecams-method
Summary: A brief discussion of Le Cam's two-point argument for minimax lower bounds.
$$
    \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    \newcommand{\EE}{\mathbf{E}}
    \newcommand{\Pp}{\mathcal{P}}
    \newcommand{\Dd}{\mathcal{D}}
$$

Le Cam's two-point method is a fundamental tool in the establishment of minimax
lower bounds.  After identifying suitable subsets of the parameter space, the
key to the argument is a lower bound of the maximum by an average.  Le Cam's
method can be found at the heart of other minimax lower bound arguments, such
as [Assouad's Lemma]({filename}assouads-lemma.md).  For this discussion, we
will follow the exposition from Bin Yu's Chapter 29 of [Festschrift for
Le Cam](http://link.springer.com/book/10.1007%2F978-1-4612-1880-7).

## Statement
Suppose $\Pp$ is a family of probability measures, with each probability
measure $P$ parameterized by a $\theta(P)$ which lives in a pseudo-metric
space $(\Dd, d)$.  

**Remark.**  Yu notes that we use pseudo-metric spaces as requiring $d(\theta,
\theta') \Leftrightarrow \theta = \theta'$ would be troublesome.

We define $\hat\theta = \hat\theta(X)$ to be the estimator of $\theta(P)$
based on $X$ drawn from a distribution $P$.  Finally, we denote the convex
hull of $\Pp$ by $co(\Pp)$.

**Lemma.** _Let $\hat\theta$ be an estimator of $\theta(P)$ on $\Pp$ taking
values in a metric space $(\Dd, d)$.  Suppose that there are subsets $\Dd_1$
and $\Dd_2$ of $\Dd$ that are $2\delta$-separated; that is, $d(s_1, s_2)
\geq 2\delta$ for all $s_1\in \Dd_1, s_2\in\Dd_2$.  Suppose further that
$\Pp_1, \Pp_2$ are subsets of $\Pp$, such that $\theta(P)\in\Dd_1$ for
$P \in \Pp_1$< and $\theta(P)\in\Dd_2$ for $P\in\Pp_2$.  Then
$$
\sup_{P\in\Pp}\EE_Pd(\hat\theta, \theta(P)) \geq
    \delta \cdot \sup_{P_1\in co(\Pp_1),\\ P_2\in co(\Pp_2)}
    \norm{P_1\wedge P_2}
$$
where $\norm{P_1\wedge P_2}$ is the total variation affinity._

_Proof._  We begin by observing that we may bound the supremem from below
by an average over an arbitrary pair of $P_1 \in \Pp_1$ and $P_2 \in \Pp_2$:
\begin{align*}
\sup_{P\in\Pp}\EE_Pd(\hat\theta, \theta(P))
    &\geq\frac{1}{2}\left[
    \EE_{P_1} d(\hat\theta, \theta(P_1))
    +
    \EE_{P_2} d(\hat\theta, \theta(P_2))
    \right]
\end{align*}
We are there able to lower bound the distances to $\theta(P_1) \in \Dd_1$ and
$\theta(P_2) \in \Dd_2$ by the distances to arbitrary members of $\Dd_1$ and
$\Dd_2$ -- say, the elements in each subset closest to $\hat\theta$.
\begin{align*}
    \frac{1}{2}\left[
    \EE_{P_1} d(\hat\theta, \theta(P_1))
    +
    \EE_{P_2} d(\hat\theta, \theta(P_2))
    \right]
    &\geq
    \frac{1}{2}\left[
    \EE_{P_1} d(\hat\theta, \Dd_1)
    +
    \EE_{P_2} d(\hat\theta, \Dd_2)
    \right]
\end{align*}
At this point, we may note that each expectation is an integration with respect
to a probability measure.  Therefore, we may interpret each one as a weighted
sum, and lower bound the sum of the two by integrating with respect to the
pointwise minimum of the two measures.
\begin{align*}
    \frac{1}{2}\left[
    \EE_{P_1} d(\hat\theta, \Dd_1)
    +
    \EE_{P_2} d(\hat\theta, \Dd_2)
    \right]
    &\geq
    \frac{1}{2}\left[
    \int \left(d(\hat\theta, \Dd_1) + d(\hat\theta, \Dd_2)\right)
        \min\{p_1, p_2\}d\mu
    \right]
\end{align*}
Applying the triangle inequality completes the proof:
\begin{align*}
    \frac{1}{2}\left[
    \int \left(d(\hat\theta, \Dd_1) + d(\hat\theta, \Dd_2)\right)
        \min\{p_1, p_2\}d\mu
    \right]
    &\geq
    \frac{1}{2}\left[
    \int 
    d(\Dd_1, \Dd_2)
        \min\{p_1, p_2\}d\mu    
    \right] \\
    &\geq
    \frac{1}{2}\left[
    \int 
    2\delta
        \min\{p_1, p_2\}d\mu
    \right] \\
    &=
    \delta
    \int 
        \min\{p_1, p_2\}d\mu    \\
    &=  \delta \norm{P \wedge Q}
\end{align*}

## Discussion
Le Cam's method gives a simple, intuitive way to give a lower bound on the
maximum risk associated with estimating a parameter by considering two
subsets of the parameter space with a suitable degree of separation.  The
true challenge is identifying these subsets of the parameter space in a way
that gives a useful bound.

Better rates may be achieved by using Assouad's Lemma, which applies Le Cam's
method in a multiple comparisons setting.  Bin Yu also remarks that "Le Cam's
method often gives the optimal rate when a real function is estimated...
[whereas Assouad and Fano] seem to be effective when the whole unknown function
is being estimated.
