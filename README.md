# nn_c103_qcotd_swave_only

This code performs the `q cot d(q)` analysis of the energy levels for the publication [arXiv:2009.11825](https://arxiv.org/abs/2009.11825).

Our web-architect started an unexpected (brief) leave of absence, so this readme is currently very crude and light-weight, as well as our accompanying webpage,  [https://laphnn.github.io/nn_c103_qcotd_swave_only/](https://laphnn.github.io/nn_c103_qcotd_swave_only/).  This will be improved soon.

## Installation
This code assumes standard scientific Python libraries are installed.  In addition, we utilize four non-standard libraries

- [gvar](https://github.com/gplepage/gvar) for handling correlated, gaussian distributions (P. Lepage)

- [lsqfit](https://github.com/gplepage/lsqfit) for performing numerical minimization (P. Lepage)

- [Two Hadrons In Box](https://github.com/ebatz/TwoHadronsInBox) `C++` code to perform analysis of two hadrons in a finite volume, accompanying the publication by C. Morningstar et al, [Nucl.Phys. B 924 (2017) 477-507](10.1016/j.nuclphysb.2017.09.014)[[arXiv:1707.05817](https://arxiv.org/abs/1707.05817)]

- [pythib](https://github.com/ebatz/pythib) A Python interface to TwoHadronsInBox

The first two can be installed with `pip`.  The latter two software suites include instructions for installation.  In order to run `qcotd_inverse_qsq.py` for example, for now, you must edit the location of the `BMat` library at the top of the script

```
BMat_path='/Users/walkloud/work/research/c51/x_files/code/pythib'
```


## Usage

We will add more analysis code here soon.  In the meantime, the phase shift analysis can be performed by running
```
python qcotd_inverse_qsq.py
```
We will add the `spline/gradient` method as well as a more sophisticated `spectrum` analysis soon.
