Spin Transport

### Model

http://ricopic.one/spin-transport

* Numerical Simulations have been done before
* Previous simulations didn't include all the phenomena
* Spins were flipped manually
* Integrated bloch-torri equations into model
* The model is a PDE

### Fenics

It is easiest to run this in a docker container. A command is now available, `fenicsproject run`, to start the docker container. The current directory is then available in the container.

### Viewing Results

The program ParaView is used to view the results. Make sure to click the `apply` button to show the plot.

The results can also be rendered using a python script, see the ParaView Guide,

https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.3&type=data&os=all&downloadFile=ParaViewGuide-5.3.0.pdf

### Multithreading

Fenics can be multithreaded using MPI. This is easy as running `mpirun -n 4 python ft01_poisson.py` where `n` is the number of threads to use and `ft01_poisson.py` is the python script to use. It appears that in most cases no extra work is needed to allow multithreading. MPI may also be used to spawn instances accross multiple networked machines.

### The System

$$\int_0^L\partial_{\bar{r}}^2\rho_i\,v_i\,d\bar{r}=-\int_0^L\partial_{\bar{r}}\rho_i\,d\bar{r}\,v_i\,d\bar{r}+v_i(L)\,\partial_{\bar{r}}\rho_i(L)-v_i(0)\,\partial_{\bar{r}}\rho_i(0)$$

### Finite Difference

https://en.wikipedia.org/wiki/Finite_difference

$$f'(x)=\frac{f(x)-f(x-h)}{h}$$

$$f''(x)=\frac{\sum_{i=0}^2(-1)^i\left(\begin{matrix}2\\i\end{matrix}\right)f(x-ih)}{h^2}$$

$$f''(x)=\frac{f(x)-2f(x-h)+f(x-2h)}{h^2}$$

### The System

$$\int_0^L\frac{\rho_i(\bar{t},\bar{r})-2\rho_i(\bar{t}-d\bar{t},\bar{r})+\rho_i(\bar{t}-2d\bar{t},\bar{r})}{d\bar{t}^2}v_i\,d\bar{r}=-\int_o^L\frac{\rho_i(\bar{t},\bar{r})-\rho_i(\bar{t}-d\bar{t},\bar{r})}{d\bar{t}}d\bar{r}\,v_i\,d\bar{r}+v_i(L)\,\partial_{\bar{r}}\rho_i(L)-v_i(0)\,\partial_{\bar{r}}\rho_i(0)$$
