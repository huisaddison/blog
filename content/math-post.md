Title: Trying Out MathJax
Date: 2017-01-25
Category: Pelican
Tags: pelican, mathjax, math
Slug: trying-out-mathjax
Summary: A post to try out MathJax

Pelican has the ability to render mathematics using
[MathJax](https://www.mathjax.org/) by way of the [Math
Render](https://github.com/getpelican/pelican-plugins/tree/master/render_math) 
plugin.

Here, we show the Cauchy-Schwarz Inequality:
$$
\langle f, g\rangle^2 \leq \langle f, f\rangle \cdot \langle g, g\rangle
$$

And here is an align environment:
\begin{align*}
w_{ij} = \begin{cases}
    1   &\text{ for } |i - j| \leq \frac{k}{2}   \\
    2 - 2\frac{|i-j|}{k} &\text{ for } \frac{k}{2} \leq |i - j| \leq k  \\
    0   &\text{ otherwise }
\end{cases}
\end{align*}

