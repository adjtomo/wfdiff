# Notes on Documentation Building

**This is only the readme for how to build the documentation. For a rendered
 view of the documentation go [here](http://krischer.github.io/veritas/).**

The complete documentation is contained in the `index.ipynb` IPython notebook.
Always edit that; the `index.rst` file will be overwritten by the sphinx
makefile.

Please don't commit the notebook with filled cells! Running `make html` will
also strip the output from the notebook. Make sure to run it at least once
before commiting changes to the notebook!


To actually build the documentation run

```bash
make html
```

This will execute the notebook and convert it to an rst file which will then
be used by sphinx to build the documentation. It will also strip the output
and some other noise from the notebook.

Install everything needed to build the documentation (in addition to veritas)
with

```bash
$ pip install sphinx sphinx-readable-theme ipython runipy
```

## Logo

The veritas logo is available in two versions:

* `logo.svg` SVG format suitable for the web.
* `logo_ref.svg` Inkscape specific svg suitable for editing in Inkscape.
