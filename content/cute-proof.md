Title: A Cute Proof (Product of Kernels is a Kernel)
Date: 2017-02-05
Modified: 2017-02-05
Category: Statistics
Tags: stats, learning-theory, rkhs, kernels
Slug: cute-proof-about-kernels 
Summary: A cute proof about kernels.

## Background
One of the most fundamental  concepts of statistical learning theory is that
of the Reproducing Kernel Hilbert Space.  Recall that a kernel function on an
RKHS is such that $K(x, y) = \langle K_x, K_y \rangle_\mathcal{H}$.  Rather
than evaluating an inner product in the Hilbert space, we may simply evaluate
the kernel function.

More generally, we may define a positive semidefinite kernel $K:\mathcal{X}
\times \mathcal{X} \mapsto \mathbf{R}$ such that
$M = \left(K(x_i, x_j)\right)_{ij}$ is positive semidefinite for all finite
collections $\{x_i\}_{i=1}^n$.

To make this result more useful, it would be nice te be able to have a
consistent way of discovering or constructing kernels.  One such way is to
take the product of two existing kernels.

## Statement & Proof
**Claim.** _Suppose $K_1, K_2: \mathcal{X}\times\mathcal{X} \mapsto \mathbf{R}$
are kernels.  Then
$$
K(x, y) = K_1(x, y) K_2(x, y)
$$
is a kernel._

_Proof._  Because $K_1, K_2$ are kernels, given an arbitrary finite collection
of $\{x_1\}_{i=1}^n$, the matrix defined by evaluating the kernel function on
all pairs of the collection satisfies the following:
$$
\mathbf{R}^{n\times n} \ni K_j = \left(K_j(x_i, x_k)\right)_{ik} \succeq 0
$$
for $j \in \{1, 2\}$.  Let $\mathbf{R}^n \ni X_1 \sim \mathcal{N}(0, K_1),
X_2 \sim \mathcal{N}(0, K_2)$ be independent.  Define $Y \triangleq X_1
\circ X_2$ the elementwise product of $X_1$ and $X_2$.  Then:
\begin{align*}
\mathrm{cov}(Y_i, Y_j) 
    &=  \mathbf{E} Y_i Y_j \\
    &=  \mathbf{E} X_{1i}X_{2i}X_{1j}X_{2j} \\
    &=  \mathbf{E} X_{1i}X_{1j}\mathbf{E}X_{2i}X_{2j} \\
    &=  \mathrm{cov}(X_{1i}, X_{1j})\mathrm{cov}(X_{2i}, X_{2j})  \\
    &=  K_1(x_i, x_j)K_2(x_i, x_j)      \\
    &=  K(x_i, x_j)
\end{align*}
and we know covariance matrices in general are positive semidefinite.

## Discussion
I'm particularly fond of this proof because it elegantly connects two facts:

*   Kernel functions evaluated on finite sets yield positive semidefinite
    matrices by definition.
*   Covariance matrices are by definition positive semidefinite.

We begin with two kernel functions, and show that the product of their outputs
may be considered valid output for an other kernel; but along the way, we treat
the matrices created by evaluating the two known kernel functions on finite
collections as covariance matrices, and construct a new covariance matrix that
coincides with the desired matrix associated with the function we sought to
prove is a kernel.  With the multivariate normal distribution as our aid, we
step in to and back out of the world of probability to prove a result in the
field of learning theory.
