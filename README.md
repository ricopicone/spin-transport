# spin-transport

This repository contains the (developing) open-source code for simulating bulk spin transport&mdash;diffusion and separation&mdash;in solid media. Multi-spin-species and magnetic resonance simulations are in development.

This is a [Python](https://www.python.org/) and [FEniCS](https://fenicsproject.org/) project. FEniCS is used to numerically solve the spin transport governing partial differential equations.

End users of this project write Python code to interface with FEniCS.

## Installation

One must first have a working installation of FEniCS.
This README assumes the use of [Docker](https://www.docker.com/) for installation, which is documented [here](http://fenics.readthedocs.io/projects/containers/en/latest/).

Next, [clone](https://help.github.com/articles/cloning-a-repository/) this repository to the host machine.

Finally, for a smooth workflow, from your working directory, create a link in your path to the `fenics` executable bash script (assuming *nix operating system).

```bash
ln fenics /usr/local/bin
```

## Workflow

## Acknowledgement

This work is supported by a grant from the Army Research Office, Materials Science Division under grant proposal _Nanoscale Spin Hyperpolarization and Imaging_ <!-- TODO: grant number -->
with PI [John Marohn, PhD](http://marohn.chem.cornell.edu/).

## Contributors

This project stems from a collaboration among three institutions:

- [Cornell University](http://www.cornell.edu/), 
- [Saint Martin's University](https://www.stmartin.edu/), and the
- [University of Washington](http://www.washington.edu/).

The lead contributor to this project is [Rico Picone, PhD](http://ricopic.one) of Saint Martin's University, co-PI on the ARO grant.
Other contributors include [John Marohn, PhD](http://marohn.chem.cornell.edu/) (Cornell, PI), John A. Sidles, PhD (Washington), Joseph L. Garbini, PhD (Washington), and Corinne Isaac (Cornell).