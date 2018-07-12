## spin-transport README

This repository contains the (developing) open-source code for simulating bulk spin transport&mdash;diffusion and separation&mdash;in solid media. Multi-spin-species and magnetic resonance simulations are in development.

This is a [Python](https://www.python.org/) and [FEniCS](https://fenicsproject.org/) project. FEniCS is used to numerically solve the spin transport governing partial differential equations.

End users of this project write Python code to interface with FEniCS.

### Installation

One must first have a working installation of FEniCS.
This README assumes the use of [Docker](https://www.docker.com/) for installation, which is documented [here](http://fenics.readthedocs.io/projects/containers/en/latest/).

Install `fenicsproject` command with the following command.

```console
curl -s https://get.fenicsproject.org | bash
```

Then [clone](https://help.github.com/articles/cloning-a-repository/) this repository to the host machine.

### Workflow

The FEniCS docs have a section on [workflow](http://fenics.readthedocs.io/projects/containers/en/latest/work_flows.html).
There are many ways to instantiate these good practices, but we recommend using the following command, which opens a Docker container, executes the `<python script>`, and saves the results to a local file (outside the Docker container).

```console
fenicsproject run stable python3 <python script> <options>
```

<!-- if you're using a \*nix system, the following may be the easiest.

With the cloned `spin-transport` repository as your working directory, create a link in your path to **spin-transport**'s `fenics` executable bash script.


```shell
ln fenics /usr/local/bin
``` -->

<!-- Now a FEniCS Python script `foo.py` can be started with the command `fenics foo.py` **from the host** instead of manually starting it from a Docker container.
This has several advantages, including that there is no need to move scripts into the container and that the complicated syntax need not be remembered. -->

### Testing the installation

To verify that everything is installed correctly, run the Poisson equation demo `ft01_poisson.py` ([source](https://fenicsproject.org/pub/tutorial/html/._ftut1004.html)) in your container.

If you installed the `fenics` bash script per the instructions above, you can use the following command (working directory: `spin-transport`).


```shell
$ fenics ft01_poisson.py
```

If everything is working fine, the output should look something like the following.


```shell
$ fenics ft01_poisson.py
Calling DOLFIN just-in-time (JIT) compiler, this may take some time.
--- Instant: compiling ---
Calling FFC just-in-time (JIT) compiler, this may take some time.
Calling FFC just-in-time (JIT) compiler, this may take some time.
Calling FFC just-in-time (JIT) compiler, this may take some time.
Solving linear variational problem.
*** Warning: Degree of exact solution may be inadequate for accurate result in errornorm.
Calling FFC just-in-time (JIT) compiler, this may take some time.
Calling FFC just-in-time (JIT) compiler, this may take some time.
  Ignoring precision in integral metadata compiled using quadrature representation. Not implemented.
Calling FFC just-in-time (JIT) compiler, this may take some time.
Calling FFC just-in-time (JIT) compiler, this may take some time.
Calling FFC just-in-time (JIT) compiler, this may take some time.
error_L2  = 0.00823509807335
error_max = 1.33226762955e-15
```

The directory `spin-transport/poisson` should have been created and should contain two files: `solution.pvd` and `solution000000.vtu`.
These files contain the solution data.

### Running a spin transport simulation

The `spin_transport_01.py` script can be executed with the following command.

```console
fenicsproject run stable python3 spin_transport_01.py -p -b
```

The default steady state IC is derived from a system with no Bloch relaxation or pulse dynamics, which are disabled with the `-p` and `-b` command line options. The number of cells in the mesh can be changed using the `-n` option followed by the number of cells to use. The full list of possible commands can be viewed by using the `-h` option.

### Plotting the results of a spin transport simulation

Use the `plot.py` script, as follows.

```console
python3 plot.py
```

Here the `-p` option can be used enable plotting with subplots. The `-i` option is used to set the duration to wait between showing the next timestep in ms. The `--start` and `--end` options are used to specify the index of the first and last timesteps to show. The rest of the options can be viewed using the `-h` option.

### Acknowledgement

This work is supported by a grant from the Army Research Office, Materials Science Division under grant proposal **Nanoscale Spin Hyperpolarization and Imaging**
with PI [John Marohn, PhD](http://marohn.chem.cornell.edu/).

### Contributors

This project stems from a collaboration among three institutions:

* [Cornell University](http://www.cornell.edu/),
* [Saint Martin's University](https://www.stmartin.edu/), and the
* [University of Washington](http://www.washington.edu/).

The lead contributor to this project is [Rico Picone, PhD](http://ricopic.one) of Saint Martin's University, co-PI on the ARO grant.
Other contributors include [John Marohn, PhD](http://marohn.chem.cornell.edu/) (Cornell, PI), John A. Sidles, PhD (Washington), Joseph L. Garbini, PhD (Washington), and Cameron Devine (Washington).
