Title: Kernel Ridge Regression
Date: 2017-02-05
Modified: 2017-02-05
Category: Statistics
Tags: stats, learning-theory, rkhs, kernels, ridge-regression
Slug: kernel-ridge-regression
Summary: A short discussion of the subtle power of kernel ridge regression.

*   The representer theorem tells us that we may express the solution to an
    optimization problem with any loss function and even a squared norm penalty
    solely in the span of the functionals defined by evaluating the kernel
    function on the feature vectors.
*   With a bit of linear algebra trickery, we are able to express the solution
    for each $X_i$ as the $i$th row of $K$ re-weighted by $\alpha$.
    Essentially, we are able to learn weights directly on the Gram matrix!
*   In doing so, we took a potentially infinite dimensional nonparametric
    problem (consider the case of Sobolev space) and turned it into a finite
    dimensional parametric problem!
