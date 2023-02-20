## Installation instructions
Assuming you have already installed R on your system (and also R studio, preferrably)
0. Clone this repo and "cd" inside it.
1. Make sure Ipopt is installed: for a tutorial go [here](https://coin-or.github.io/Ipopt/INSTALL.html). If you are using **arch linux**, you can simply run `yay -S coin-or-ipopt`
2. Make sure CppAD is installed, instructions [here](https://coin-or.github.io/CppAD/doc/install.htm). 
<br>
I recommend you perform the  `git clone ` command in the working directory of this project, since you might have to modify the -I directive in the `src/Makevars` file otherwise.
<br>
**Important: remember to add the option `-D include_ipopt=true` at step 2**, see [here]{https://www.coin-or.org/CppAD/Doc/ipopt_prefix.htm}, else you will get a missing header file..

3. Compile and install the quadrature library provided by **pacs-examples** (see the repo for requirements)
  * Go to `$(PACS_ROOT)/src/QuadratuleRule 
  * Follow the instructions on the README.md file to compile and create the .so of the baseVersion
  **Nota bene**: you will need to have the [OpenBlas lib](https://github.com/xianyi/OpenBLAS/wiki/Installation-Guide) installed beforehand.
  
4. Go to `.src/Makevars` and set the first variables that are signalled to be edited
* CPPAD_LIB_DIR with the directory where the .so file of cppad was installed; CPPAD_INC_DIR for its header files
* IPOPT_LIB_DIR for the .so file of ipopt; IPOPT_INC_DIR for its headers.
* PACS_ROOT the path we saw at the pacs course, i.e. if you clone the repo here[https://github.com/pacs-course/pacs-examples], its path so that you are inside the `Examples` folder.
* mkOpenBlasLib the path to the library OpenBlas (however this is usually on the directories where programs are regularly installed and you should not need to add it.)
<br> 
If you do not like editing the `Makevars` file,, you can do the following:
    * (optional) Set the environment variable R_HOME, which may be useful for you later in your life.
    * Run `R -e "R.home()"` to find the value of R_HOME
    * Go to **R_HOME/etc/Renviron** and set the variables there.
5. Go into the `inst/` directory and clone [this repo]{https://github.com/gcant/dirichlet-cpp/tree/4023dc18599a621c7452caf8be16494dcc00ee2a}, and go back to the directory of FdPot (inside FdPot, outside inst/)
6. Now, you can run `make install` to install the package

7. You are all set, open `R` and to load the library> `library(FdPot)`

8. Run `?FdPot` for help

### Other features.
After installation, run: `make docs` to obtain the documentation, it will be created inside the `.src/doc` folder.
  