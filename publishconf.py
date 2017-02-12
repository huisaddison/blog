#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals

# This file is only used if you use `make publish` or
# explicitly specify it as your config file.

import os
import sys
sys.path.append(os.curdir)
from pelicanconf import *

HOMEURL = 'http://huisaddison.com'
HOMENAME = 'huisaddison'

SITEURL = 'http://huisaddison.com/blog'
SITENAME = 'huisaddison/blog'


RELATIVE_URLS = False

DELETE_OUTPUT_DIRECTORY = True

PLUGIN_PATHS = ['../pelican-plugins']
PLUGINS = ['render_math']

THEMES = '../blog-theme'


# Following items are often useful when publishing

#DISQUS_SITENAME = ""
#GOOGLE_ANALYTICS = ""
