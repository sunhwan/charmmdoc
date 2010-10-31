import re
from docutils import nodes, utils

_s = re.compile('\s+')

def charmm_keyword_role(name, rawtext, text, lineno, inliner, option={}, content=[]):
    keywords = _s.split(text)
    text = []
    for kwd in keywords:
        if len(kwd) > 4:
            text.append(kwd[:4].upper() + kwd[4:].lower())
        else:
            text.append(kwd.upper())
    node = nodes.literal(rawtext, ' '.join(text))
    return [node], []

def setup(app):
    app.add_role('chm', charmm_keyword_role)
