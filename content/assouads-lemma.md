Title: Assouad's Lemma
Date: 2017-02-01
Modified: 2017-02-01
Category: Statistics
Tags: math, stats, minimax, risk, assouad
Slug: assouads-lemma
Summary: A brief discussion of Assouad's Lemma, its uses, and its relationship to Le Cam's Method.
$$
    \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    \newcommand{\EE}{\mathbf{E}}
$$
Assouad's Lemma is a common tool for proving lower bounds on the maximum risk
associated with an estimation problem, and its core it may be interpreted as
a multiple comparisons version of Le Cam's Method.  This post will first
review the version discussed by Bin Yu in her chapter of [Festschrift
for Lucien Le Cam](http://www.springer.com/us/book/9780387949529).  After that,
we will discuss a variation found in Aad van der Vaart's [Asymptotic
Statistics](
https://books.google.com/books/about/Asymptotic_Statistics.html?id=UEuQEM5RjWgC
).  We conclude with a brief discussion of its use.

## Assouad à la Bin Yu
### Preliminaries

*   Assume that $\mathcal{P}$ is a family of probability measures and $\theta(P)$
    is the parameter of interest with values in the pseudo-metric space
    $(\mathcal{D}, d)$.
*   Let $\hat\theta = \hat\theta(X)$ be an estimator based on an $X$ drawn
    from a distribution $P$.

### Statement
Bin Yu gives the following statement of Assouad's Lemma in Chapter 29:

**Lemma 2 (Assouad's Lemma)**  _Let $m\geq 1$ be an integer and let
$\mathcal{F}_m = \{P_\tau: \tau \in \{-1, +1\}^m\}$ contain $2^m$ probability
measures.  Write $\tau \sim \tau'$ if $\tau$ and $\tau'$ differ in only one
coordinate, and write $\tau \sim_j \tau'$ when that coordinate is the jth.
Suppose that there are $m$ pseudo-distances on $\mathcal{D}$ such that
for any $x, y \in \mathcal{D}$
$$
    d(x, y) = \sum_{j=1}^m d_j(x, y)
$$
and further, that if $\tau\sim_j\tau'$,
$$
    d_j(\theta(P_\tau), \theta(P_{\tau'})) \geq \alpha_m
$$
Then_
$$
    \max_{P_\tau\in\mathcal{F}_m}
    \mathbf{E}_\tau d(\hat\theta, \theta(P_\tau)) \geq
    m\cdot \frac{\alpha_m}{2}\min\{\norm{P_\tau\wedge P_{\tau'}}:
    \tau\sim\tau'\}
$$

### Proof
We now discuss the proof provided.  For each tuple $\tau = 
(\tau_1, \dotsc, \tau_m)$, let $\tau^j$ denote the $m$-tuple
that differs only at the jth index.  Then $d(\theta(\tau),
\theta(P_{\tau^j}))\geq \alpha_m$ by assumption.

The key idea in this proof that we may bound the maximum risk over $\tau$
from below by the average over the risk for each $\tau$, and then
apply Le Cam's Method many times to get our desired bound.

But first, we use the assumption that we may decouple our pseudo-distance
metric across the members of $\tau$:
\begin{align*}
    \max_\tau \EE_\tau d(\theta(P_\tau), \hat\theta)
        &=  \max_\tau\sum_{j=1}^m \EE_{\tau_j} d(\theta(P_\tau), \hat\theta)
\end{align*}
after which we use the average to lower bound and rearrange the summations:
\begin{align*}
    \max_\tau \EE_\tau d(\theta(P_\tau), \hat\theta)
        &=  \max_\tau\sum_{j=1}^m \EE_{\tau_j} d(\theta(P_\tau), \hat\theta)\\
        &\geq 2^{-m}\sum_\tau \sum_{j=1}^m
            \EE_\tau d_j(\theta(P_\tau), \hat\theta)    \\
        &= \sum_{j=1}^m2^{-m}\sum_\tau 
            \EE_\tau d_j(\theta(P_\tau), \hat\theta)
\end{align*}
afterwards, Yu cleverly rearranges the terms so that each $\tau$ is matched
up with a $\tau^j$ by strategically adding a copy of each term and then
dividing through by a half:
\begin{align*}
    \sum_{j=1}^m2^{-m}&\sum_\tau 
        \EE_\tau d_j(\theta(P_\tau), \hat\theta)    \\
    &=  
    \sum_{j=1}^m2^{-m}\sum_\tau 
        \frac{1}{2}
        \left(
        \EE_\tau d_j(\theta(P_\tau), \hat\theta)    +
        \EE_{\tau^j} d_j(\theta(P_{\tau^j}), \hat\theta)
        \right)
\end{align*}

For each $\tau$ and $j$, consider the pair of associated hypotheses $P_\tau$
and $P_{\tau^j}$.  By using Le Cam's Method, we may bound the average
estimation error for each of these pairs from below:
\begin{align*}
    \max_\tau \EE_\tau d(\theta(P_\tau), \hat\theta)
        &\geq \sum_{j=1}^m 2^{-m} \sum_\tau
            \frac{\alpha_m}{2} \norm{P_\tau \wedge P_{\tau^j}}    \\
        &\geq m
            \frac{\alpha_m}{2}
             \min\{\norm{P_\tau \wedge P_{\tau^j}}: \tau\sim\tau'\}
\end{align*}
giving our desired bound.

## Assouad à la Aad van der Vaart
The [paper](https://arxiv.org/abs/1010.3866) that I've been reading lately by
Zhou et al uses a slightly different version of Assouad's Lemma, elaborated by
van der Vaart in his book Asymptotic Statistics.  My post about lower bounds
discussed in that paper (which takes advantage of Assouad's Lemma) may be found
[here]({filename}discuss-covmat-optimal-convergence-pt2.md).

The version discussed by van der Vaart relaxes the assumption that the distance
metric decouple across $j$.

### Preliminaries
Suppose we have a parameter set $\Theta = \{0, 1\}^r$ and are estimating
an arbitrary quantity $\psi(\theta)$, belonging to a metric space with
metric $d$.

### Statement
**24.3  Lemma (Assouad).**  _For any estimator $T$ based on an observation
in the experiment ($P_\theta: \theta\in\{0, 1\}^r$), and any $p > 0$,_
$$
\max_\theta 2^p \EE_\theta d^p(T, \psi(\theta))
\geq
\min_{H(\theta, \theta') \geq 1} \frac{
    d^p(\psi(\theta), \psi(\theta'))
}{
    H(\theta, \theta')
}
\frac{r}{2}
\min_{H(\theta, \theta') = 1} \norm{P_\theta \wedge P_{\theta'}}
$$

We see that the first term in the lower bound, which gives the minimum
distance between two parameters indexed by the vertices of a hypercube,
is normalized by the Hamming distance between the vertex indices, whereas
in Yu's version of Assouad's Lemma the assumption that we could decompose
the distance metric allowed us to only look at parameters corresponding
to neighboring vertices.

### Proof
We define an estimator $S$ as follows:
$$
S \in\arg\min_{\{0, 1\}^r} d(T, \psi(S))
$$
For any $\theta$, we have the lower bound:
\begin{align*}
d(\psi(S), \psi(\theta)
    &\leq d(\psi(S), T) + d(\psi(\theta), T)
    &&\text{(Triangle Inequality.)} \\
    &\leq 2d(\psi(\theta), T)
    &&\text{(By definition of $S$)}
\end{align*}

Pick $\alpha$ such that for all $\theta \neq \theta'$:
$$
d^p(\psi(\theta), \psi(\theta')) \geq \alpha H(\theta, \theta')
$$
Equivalently,
$$
\alpha = \min_{\theta \neq \theta'}
\frac{d^p(\psi(\theta), \psi(\theta'))}{H(\theta, \theta')}
$$

Therefore, we may bound:
\begin{align*}
    2^p\EE_\theta d^p(T, \psi(\theta))
    &\geq 2^p \EE_\theta \frac{1}{2^p}d^p(\psi(S), \psi(\theta))    \\
    &\geq \alpha \EE_\theta H(S, \theta)
\end{align*}
Let's focus on the expectation term and bound the maximum by the average.
\begin{align*}
\max_\theta \EE_\theta H(S, \theta)
&\geq  \frac{1}{2^r}\sum_\theta \EE_\theta H(S, \theta)     \\
&      \frac{1}{2^r}\sum_\theta 
    \sum_{j=1}^r \EE_\theta|S_j-\theta_j|       \\
&      \frac{1}{2^{r-1}}
    \sum_{j=1}^r
    \sum_{\theta: \theta_j \in\{0, 1\}}
    \frac{1}{2}
    \left(
    \EE_{\theta: \theta_j = 0} |S_j-\theta_j|   +
    \EE_{\theta: \theta_j = 1} |S_j-\theta_j|       \right) \\
&=     \frac{1}{2^{r-1}}
    \sum_{j=1}^r
    \sum_{\theta: \theta_j \in\{0, 1\}}
    \frac{1}{2}
    \left(
    \int (1 - S_j) p_{\theta:\theta_j = 0}\;d\mu   +
    \int S_j p_{\theta: \theta_j = 1}\;d\mu        \right) \\
&\geq  \frac{1}{2^{r-1}}
    \sum_{j=1}^r
    \sum_{\theta: \theta_j \in\{0, 1\}}
    \frac{1}{2}
    \left(
    \int (1 - S_j) \min\{p_{\theta:\theta_j = 0},
        p_{\theta: \theta_j = 1}\}\;d\mu   +
    \int S_j \min\{p_{\theta:\theta_j = 0},
        p_{\theta: \theta_j = 1}\}\;d\mu \right) \\
&=     
    \frac{1}{2}
    \sum_{j=1}^r
    \frac{1}{2^{r-1}}
    \sum_{\theta: \theta_j \in\{0, 1\}}
    \int \min\{p_{\theta:\theta_j = 0},
        p_{\theta: \theta_j = 1}\}\;d\mu  \\
&\geq
    \frac{1}{2}
    \sum_{j=1}^r
    \frac{1}{2^{r-1}}
    \sum_{\theta: \theta_j \in\{0, 1\}}
    \norm{P_{\theta:\theta_j = 0} \wedge
        P_{\theta: \theta_j = 1}}\\
&\geq
    \frac{1}{2}
    \sum_{j=1}^r
    \min_{H(\theta, \theta') = 1}
    \norm{P_{\theta} \wedge
        P_{\theta'}}            \\
&=  \frac{r}{2}
    \min_{H(\theta, \theta') = 1}
    \norm{P_{\theta} \wedge
        P_{\theta'}}
\end{align*}

Appending the factor of $\alpha$ previously calculated, we have our desired
bound:
\begin{align*}
\max_\theta 2^p \EE_\theta d^p(T, \psi(\theta))
&\geq
\alpha
\frac{r}{2}
\min_{H(\theta, \theta') = 1} \norm{P_\theta \wedge P_{\theta'}}    \\
&\geq
\min_{H(\theta, \theta') \geq 1} \frac{
    d^p(\psi(\theta), \psi(\theta'))
}{
    H(\theta, \theta')
}
\frac{r}{2}
\min_{H(\theta, \theta') = 1} \norm{P_\theta \wedge P_{\theta'}}
\end{align*}

## Discussion
As previously observed, Assouad's Lemma may be decomposed into a set of
two-point comparisons for which Le Cam's Method may be (and is) applied.
So why bother with Assouad?

By simultaneously considering, say, $m$ two-point comparisons, we are able to
push up our lower bound by a factor of $m$ corresponding to the dimensionality
of the associated hypercube, which can be convenient or necessary to establish
the optimality of an estimator; it certainly doesn't hurt to have a tighter
bound.  As Yu remarks in her note, "it is known that Assouad's Lemma (Lemma 2)
gives very effective lower bounds for many global estimation problems." 

Of course, we must satisfy the assumptions on the distance between parameters
in our space, as well as work to construct an intelligent mapping from the
vertices of the hypercube to the parameters themselves so that the pairwise
application of Le Cam succeeds.  Regardless, it's delightful to see that Le
Cam's Method can be found at the heart of such a useful tool.
