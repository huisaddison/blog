Title: Reproducing Kernel Hilbert Spaces
Date: 2017-01-30
Category: Statistics
Tags: math, stats, support-vector-machines, kernel, hilbert-spaces
Slug: reproducing-kernel-hilbert-spaces
Summary: Discussion of the magic behind what makes support vector machines so powerful.

## Introduction

## Background

## Reproducing Kernel Hilbert Spaces

$k(\cdot, \cdot)$ is a reproducing kernel of a Hilbert space $\mathcal{H}$ if
$\forall f \in \mathcal{H}, f(x) = \langle k(x, \cdot), f(\cdot)\rangle$.

### Example: $\mathbf{R}^n$
As a simple example, let us consider the Hilbert spaces $\mathbf{R}^n$ equipped
with the inner product $\langle x, y \rangle = \sum_{i=1}^n x_i y_i$.

**Claim.** The kernel $k(x, y) = x^\top\mathbf{I}y$ is a reproducing kernel for
this Hilbert space.

*Proof.* Consider an arbitrary function $f(y) = v^\top y$.
\begin{align*}
    \langle k(x, \cdot), f(\cdot) \rangle
    &=  \langle x^\top\mathbf{I}, v \rangle    \\
    &=  \langle x, v \rangle    \\
    &=  x^\top v   \\
    &=  f(x) 
\end{align*}
<div style="text-align: right">&#8718;</div>


## Implications
The example in $\mathbf{R}^n$ is contrived, but the existence of reproducing
kernels is more powerful when the result is not so obvious.  Rather than
directly applying a function $f \in \mathcal{H}$ to $x$, we may instead apply
the kernel $k(\cdot, \cdot)$ to $x$, and then evaluate the inner product
$\langle k(x, \cdot), f(\cdot)\rangle$.


