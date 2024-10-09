#!/home/ty/anaconda3/bin/python

import argparse
import os

# --------------------------------------------------------------------------------------------------

def _get_args():

    description = 'compile a latex paper w/ bibliography or clean dir'
    cmd_parser = argparse.ArgumentParser(description=description)

    help_msg = 'compile pdf from paper w/ this prefix(s)'
    cmd_parser.add_argument('-i','--file-prefix',default=['paper'],help=help_msg,nargs='+')
    help_msg = 'clean up the directory: remove log files etc. for given prefix(s)'
    cmd_parser.add_argument('-c','--clean',action='store_const',help=help_msg,const=True,
                            default=False)
    help_msg = 'do fast compile w/o running bibtex and pdflatex again to add refs'
    cmd_parser.add_argument('-f','--fast',action='store_const',help=help_msg,const=True,
                            default=False)
    help_msg = 'use bibtex instead of biber to compile bibliography'
    cmd_parser.add_argument('-b','--bibtex',action='store_const',help=help_msg,const=True,
                            default=False)

    cmd_args = cmd_parser.parse_args()
    prefixes = cmd_args.file_prefix
    clean = cmd_args.clean
    fast = cmd_args.fast
    bibtex = cmd_args.bibtex

    return prefixes, clean, fast, bibtex

# --------------------------------------------------------------------------------------------------

def _compile(prefix,fast,bibtex):

    if bibtex:
        _bib = 'bibtex'
    else:
        _bib = 'biber'
    
    if fast:
        _cmd = f'pdflatex {prefix}.tex && rifle {prefix}.pdf'
    else:
        _cmd = f'pdflatex {prefix}.tex && {_bib} {prefix} && pdflatex {prefix}.tex ' \
               f' && pdflatex {prefix}.tex && rifle {prefix}.pdf'

    print(_cmd)
    os.system(_cmd)

# --------------------------------------------------------------------------------------------------

def _clean(prefix):

    _cmd = f'rm {prefix}.aux  {prefix}.bcf  {prefix}.log  {prefix}.out {prefix}.bbl  ' \
           f'{prefix}.blg {prefix}.run.xml'
    print(_cmd)
    os.system(_cmd)

# --------------------------------------------------------------------------------------------------

prefixes, clean, fast, bibtex = _get_args()

for prefix in prefixes:
    if clean:
        _clean(prefix)
    else:
        _compile(prefix,fast,bibtex)
    
    


