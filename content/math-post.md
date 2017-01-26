Title: Trying Out MathJax
Date: 2017-01-25
Category: Mathematics
Tags: pelican, mathjax, math
Slug: trying-out-mathjax
Summary: A post to try out MathJax

Pelican has the ability to render mathematics using
[MathJax](https://www.mathjax.org/) by way of the [Math
Render](https://github.com/getpelican/pelican-plugins/tree/master/render_math) 
plugin.

Here, we show the Cauchy-Schwarz Inequality:
$$
\langle f, g\rangle^2 = \langle f, f\rangle \cdot \langle g, g\rangle
$$

And here is an align environment:
\begin{align*}
    N^{m}
    &=  \Vert M_l^{m} - \mathbf{E}M_l^{m}\Vert  \\
    &=  \Vert \breve\Sigma - \mathbf{E}\breve\Sigma\Vert  \\
\end{align*}
