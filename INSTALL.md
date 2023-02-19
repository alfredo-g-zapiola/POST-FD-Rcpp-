## Installation instructions
Assuming you have already installed R on your system:

1. Make sure Ipopt is installed: for a tutorial go [here](https://coin-or.github.io/Ipopt/INSTALL.html), 
2. Make sure CppAD is installed, instructions [here](https://coin-or.github.io/CppAD/doc/install.htm). **Important: remember to add the option `-D include_ipopt=true` at step 2**.

3. Compile and install the quadrature library provided by **pacs-examples** (see the repo for requirements)
  * Go to `$(PACS_ROOT)/src/QuadratuleRule 
  * Follow the instructions on the README.md file to compile and create the .so of the baseVersion
  **Nota bene**: you will need to inst
  
4. Go to `.src/Makevars` and set the following variables:
* CPPAD_LIB_DIR with the directory where the .so file of cppad was installed.
* PACS_ROOT the path we saw at the pacs course, i.e. if you clone the repo here[https://github.com/pacs-course/pacs-examples], its path so that you are inside the `Examples` folder.
* mkOpenBlasLib the path to the library OpenBlas (however this is usually on the directories where programs are regularly installed and you should not need to add it.)
<br> Alternatively, you can do the following:
  a.(optional) Set the environment variable R_HOME, which may be useful for you later in your life.
  b. Run `R -e "R.home()"` to find the value of R_HOME
  c. Go to **R_HOME/etc/Renviron** and set the variables there.

  