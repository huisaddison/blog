Title: Toying with Fonts in LaTeX
Date: 2017-02-23
Modified: 2017-02-23
Category: LaTeX
Tags: latex, xelatex, fonts, fontspec
Slug: fonts-in-latex
Summary: Messing around with fonts in LaTeX

I typeset most of my documents using some program built atop TeX.  For
documents where fine control over formatting is not necessarily (e.g.,
blog posts), I use Markdown instead.  On occasion, I use LibreOffice to
type out a quick sign.

I'm particularly fond of LaTeX and its family because it does a fine job of
separating content creation from layout, and because math is rendered
beautifully (though Markdown compiled to HTML with MathJax does very well too).

Up until now, I've used pdflatex to compile my tex files, but recently I've
started using xelatex instead, in some cases.

## Fonts in LaTeX
Aside from some pre-installed fonts obtained through a distribution of LaTeX,
there are not many options for fonts when using pdflatex.  Typically this
would not be a problem, as the default font family (Computer Modern), which
was last updated in 1992[^cm-font], is quite nice.  In some cases, though,
it would be nice to be able to add some personal flair.

### fontspec
I use the moderncv package to typeset my resume, but found the default sans
serif font to be a bit bland.  I decided to replace the section headers
with [Bitter](https://fonts.google.com/specimen/Bitter) and the body text
with [Raleway](https://fonts.google.com/specimen/Raleway) by using the
fontspec package.  Defining the fonts could not be easier:
```
\usepackage{fontspec}
\newfontfamily\bitter[Path=fonts/Bitter/, Ligatures=TeX]{Bitter-Regular}
\newfontfamily\raleway[Path=fonts/Raleway/, Ligatures=TeX]{Raleway-Regular}
```
After defining the font families, I set the body font to Raleway:
```
\setsansfont[Path=fonts/Raleway/,
    BoldFont=Raleway-Bold,
    ItalicFont=Raleway-Italic,
    BoldItalicFont=Raleway-BoldItalic,
    Ligatures=TeX]{Raleway-Regular}
```
and finally, set the rest of the fonts to Bitter:
```
\renewcommand*{\namefont}{\bitter\fontsize{34}{36}\mdseries\upshape}
\renewcommand*{\titlefont}{\bitter\LARGE\mdseries\slshape}
\renewcommand*{\addressfont}{\bitter\small\mdseries\slshape}
\renewcommand*{\quotefont}{\bitter\large\slshape}
\renewcommand*{\sectionfont}{\bitter\Large\mdseries\upshape}
\renewcommand*{\subsectionfont}{\bitter\large\mdseries\upshape}
```

## XeLaTeX
The only catch is that the fontspec package is not available when using
pdflatex; the only catch is that one must use xelatex instead.  xetex
has a few advantages over pdftex[^xetex]:

* xetex assumes the input is unicode, so characters with accents and other
  marks can be inserted directly into the tex source.
* As seen in the example above, xetex is able to use fonts located on the sytem
  without any issues.


For now, I use xetex to typeset my [resume](https://github.com/huisaddison/cv)
and [slides](https://github.com/huisaddison/cpsc453-presentation); for other
documents, I stick with pdftex.

## LuaTeX
It's worth mentioning that LuaTeX is the anointed successor to pdftex, so
I will probably be using that at some point in my life.  It has ambitious
goals, like introducing Lua to TeX to allow scripting within the
document[^luatex].  I haven't tried it yet as xetex does everything I need
for now.

[^cm-font]: http://www-cs-faculty.stanford.edu/~uno/cm.html
[^xetex]: http://tex.stackexchange.com/questions/3393/
[^luatex]: http://tex.stackexchange.com/questions/36/
