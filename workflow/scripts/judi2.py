
JUDI_VALUE_SEP = "~"
JUDI_PARAM_SEP = "/"


def parse_kv_pairs(text, item_sep=",", value_sep="="):
    """Parse key-value pairs from a shell-like text."""
    from shlex import shlex
    # initialize a lexer, in POSIX mode (to properly handle escaping)
    lexer = shlex(text, posix=True)
    # set ',' as whitespace for the lexer
    # (the lexer will use this character to separate words)
    lexer.whitespace = item_sep
    # include '=' as a word character 
    # (this is done so that the lexer returns a list of key-value pairs)
    # (if your option key or value contains any unquoted special character, you will need to add it here)
    lexer.wordchars += value_sep
    # then we separate option keys and values to build the resulting dictionary
    # (maxsplit is required to make sure that '=' in value will not be a problem)
    return dict(word.split(value_sep, maxsplit=1) for word in lexer if value_sep in word)

class SafeDict(dict):
     def __missing__(self, key):
         return '{' + key + '}'

class Filespace:
    """A wrapper for pandas dataframes that provides helpers for using them as a parameter
    space in Snakemake.
    This is heavily inspired by @soumitrakp work on JUDI (https://github.com/ncbi/JUDI).
    pals2@lmem14 scripts$ python
    >>> f = Filespace("a.b")
    >>> f = Filespace("x/a.b") # error
    >>> f = Filespace("a.b", topdir="haha")
    >>> f = Filespace("a.b", topdir="haha", params=["y", "x"]).pat
    """
    def __init__(self, fname, params = None, topdir = None, name = None):
        import os, sys
        if os.path.basename(fname) != fname:
            sys.exit("fname contains directory! set topdir instead")
        if not fname:
            sys.exit("fname cannot be empty string!")
        prefix, ext = os.path.splitext(fname)
        self.fname = fname
        self.name = name or "_".join([prefix, ext[1:]])
        self.params = sorted(params or [])
        self.basedir = topdir or ""
        tmpl = "{0}"+JUDI_VALUE_SEP+"{{{0}}}"
        self.pattern = JUDI_PARAM_SEP.join(map(tmpl.format, self.params))

    def fix(self, **kwargs):
        from copy import deepcopy
        other = deepcopy(self)
        other.pattern = other.pattern.format_map(SafeDict(kwargs))
        return other


    @property
    def path(self):
        """Wildcard pattern over all columns of the underlying dataframe of the form
        column1~{column1}/column2~{column2}/***
        """
        parts = []
        if self.basedir:
            parts.append(self.basedir)
        if self.params:
            parts.append(self.name)
            parts.append(self.pattern)
        parts.append(self.fname)
        return "/".join(parts)
#
#    @property
#    def instance_patterns(self):
#        """Iterator over all instances of the parameter space (dataframe rows),
#        formatted as file patterns of the form column1~{value1}/column2~{value2}/...
#        """
#        return (
#            "/".join("{}~{}".format(name, value) for name, value in row.items())
#            for index, row in self.dataframe.iterrows()
#        )
#
#    def instance(self, wildcards):
#        """Obtain instance (dataframe row) with the given wildcard values."""
#        import pandas as pd
#
#        return {
#            name: pd.Series([value]).astype(self.dataframe.dtypes[name])
#            for name in wildcards.items()
#            if name in self.dataframe.columns
#        }
#
#    def partially_expanded_patterns(self, wildcards):
#        """Iterator over all instances of the parameter space (dataframe rows),
#        that are filtered by the wildcards values, formatted as file patterns
#        of the form column1~{z1}/column2~{z2}/... where z1 is {column1} if
#        column1 is in wildscards, otherwise value1
#        """
#        fixed = {name:value for name, value in wildcards.items()
#                 if name in self.dataframe.columns}
#        qstring = '&'.join(f"({name} == {value})"
#                           for name, value in fixed.items())
#        tmp_dataframe = (self.dataframe if not qstring else
#                         self.dataframe.query(qstring))
#        return (
#            "/".join(f"{name}~{{{name}}}" if name in fixed else
#                     f"{name}~{value}"  for name, value in row.items())
#            for index, row in tmp_dataframe.iterrows()
#        )
#
#
#    #def __getattr__(self, name):
#    #    import pandas as pd
#
#    #    ret = getattr(self.dataframe, name)
#    #    if isinstance(ret, pd.DataFrame):
#    #        return Paramspace(ret)
#    #    return ret
#
#    #def __getitem__(self, key):
#    #    import pandas as pd
#
#    #    ret = self.dataframe[key]
#    #    if isinstance(ref, pd.DataFrame):
#    #        return Paramspace(ret)
#    #    return ret
#
#
#    def expand_wildcards(self, wildcards):
#        if wildcards is None:
#            return "/".join([self.basedir,
#                             self.paramspace.wildcard_pattern,
#                             self.fname])
#        return (
#            "/".join([self.basedir, pattern, self.fname])
#            for pattern in self.paramspace.partially_expanded_patterns(wildcards)
#        )
#
#    @property
#    def expand(self):
#        return self.expand_wildcards({})
#
#    def __call__(self, wildcards=None):
#        return self.expand_wildcards(wildcards)
#
