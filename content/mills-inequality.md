Title: Mill's Inequality
Date: 2017-02-22
Modified: 2017-02-23
Category: Statistics
Tags: stats, math, concentration-bounds, gaussian, tail-bounds
Slug: mills-inequality
Summary: A useful concentration bound on Gaussian random variables.
<div style="display:none">
    $$
    \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    \newcommand{\one}{\mathbf{1}}
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

## Background
One of the results in the [paper](https://arxiv.org/abs/1309.6024) I've been
reading uses a tail bound on the max of Gaussian-distributed random variables
that I was not too familiar with, so I thought I'd present and discuss it here
to solidify my understanding of it.  

## Statement
Suppose $Z \sim \Nn(0, \sigma^2)$.  Then we have the tail bound:
$$
\PP\{|Z| > t\} \leq \sqrt\frac{2}{\pi}\frac{\sigma}{t} \exp\left\{
-\frac{t^2}{2\sigma^2}
\right\}
$$
Intuitively, this kind of tail bound is useful because we can get
exponentially-fast decay without calculating the distribution function
directly.

## Proof
The broad strokes of the proof follow [Aliyah Ahmed's response to a post
on StackExchange](http://math.stackexchange.com/questions/723041/).  We begin
by observing that density of $Z$ is symmetric about the origin, therefore:
\begin{align*}
\PP\{|Z| > t\}
&=  2 \PP\{Z > t\}
\end{align*}

We then observe that by playing with distribution functions and expectations,
we get the following upper bound:
\begin{align*}
t\cdot \PP\{Z > t\}
&= t\int_t^\infty dF(x)  \\
&\leq  \int_t^\infty x d F(x)  \\
&=  \int_t^\infty x \cdot \frac{1}{\sqrt{2\pi}\sigma}\exp\left\{
-\frac{x^2}{2\sigma^2}
\right\}    \\
&=  \frac{\sigma}{\sqrt{2\pi}}\exp\left\{
-\frac{t^2}{2\sigma^2}
\right\}
\end{align*}
in the process using sneaky way to introduce a quantity that has a nice, clean
closed-form integral.  Closer examination shows that this is in fact a tighter
version of Markov's Inequality; rather than taking $\EE X$, we take $\EE [X
\one\{X > t\}]$.  This implies that:
\begin{align*}
\PP\{Z > t\}
&=  \frac{\sigma}{t\sqrt{2\pi}}\exp\left\{
-\frac{t^2}{2\sigma^2}
\right\}    \\
\Rightarrow
\PP\{|Z| > t\}
&=  \sqrt\frac{2}{\pi} \frac{\sigma}{t}\exp\left\{
-\frac{t^2}{2\sigma^2}
\right\}
\end{align*}


## Extension to Sum of Random Variables
This result can be extended to the maximum of $m$ Gaussian random variables
by way of the union bound.  Suppose $\{Z_i\}_{i=1}^m \sim \Nn(0, \sigma^2)$.
Then the union bound implies:
$$
\PP\left\{
\max_{1\leq i\leq m} |Z_i| > t
\right\} \leq m\cdot \sqrt\frac{2}{\pi} \frac{\sigma}{t}\exp\left\{
-\frac{t^2}{2\sigma^2}
\right\}
$$
Suppose the variance of these random variables decreased with $n$, i.e.,
$\sigma^2 = \frac{1}{n}$.  This could happen if our $Z_i$ are estimators.
Then we would have the bound:
\begin{align*}
\PP\left\{
\max_{1\leq i\leq m} |Z_i| > t
\right\}
&\leq m\cdot \sqrt\frac{2}{\pi} \frac{1}{\sqrt{n}t}\exp\left\{
-\frac{nt^2}{2}
\right\}    \\
&\leq \sqrt\frac{2}{\pi} \frac{1}{\sqrt{n}t}\exp\left\{
-\frac{nt^2}{2} + \log m
\right\}
\end{align*}
Suppose we wanted to get a reasonably large probability on the right hand
side so that our bound is useful.  A trivial way to do this is to take
$t$ very small, but this upper bound would be meaningless.

What if we made $t$ arbitrarily large?  The implication in this case is not
particularly useful either.

To balance between these interests, we can choose $t$ such that:
$$
\frac{nt^2}{2} = \log m
$$
giving us a bound that adapts to $m$.
