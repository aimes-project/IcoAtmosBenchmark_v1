# What is this ? #

**IcoAtmosBenchmark** is the benchmark suite of compulational kernels from
three icosahedral atmospheric model developed and used in world wide,
**NICAM**, **DYNAMICO** and **ICON**.

Kernels are extracted to be individual executables and you are able to
run them separately.

# kernels #

In this package, kernels from **NICAM** are as belows.
See the description manual for the details of each kernel.

* `communication`
* `dyn_diffusion`
* `dyn_divdamp`
* `dyn_horiz_adv_flux`
* `dyn_horiz_adv_limiter`
* `dyn_metrics`
* `dyn_vert_adv_limiter`
* `dyn_vi_rhow_solver`

# Directory tree #

Kernel programs are stored in individual directories,
under each model's directory tree, such as
 `${TOPDIR}/kernels/NICAM/dyn_diffusion`.

Under each kernel directory, there are directories named
`data/`, `reference`, `run/` and `src/`.
They are for input/validation data files, log files for references,
execution scripts, and source files, respectively.

See below for whole directory tree.



    IcoAtmosBenchmark/
    +---- bin
    +---- doc
    |     +---- NICAM
    |     +---- DYNAMICO
    |     +---- ICON
    |
    +---- kernels
    |     +---- DYNAMICO
    |     |     +---- .....
    |     |     +---- .....
    |     +---- ICON
    |     |     +---- .....
    |     |     +---- .....
    |     +---- NICAM
    |           +---- communication
    |           |     +---- data
    |           |     +---- reference
    |           |     +---- run
    |           |     +---- src
    |           +---- dyn_diffusion
    |           |     +---- data
    |           |     +---- reference
    |           |     +---- run
    |           |     +---- src
    |           +---- dyn_divdamp
    |           |     +---- data
    |           |     +---- reference
    |           |     +---- run
    |           |     +---- src
    |           +---- dyn_horiz_adv_flux
    |           |     +---- data
    |           |     +---- reference
    |           |     +---- run
    |           |     +---- src
    |           +---- dyn_horiz_adv_limiter
    |           |     +---- data
    |           |     +---- reference
    |           |     +---- run
    |           |     +---- src
    |           +---- dyn_metrics
    |           |     +---- data
    |           |     +---- reference
    |           |     +---- run
    |           |     +---- src
    |           +---- dyn_vert_adv_limiter
    |           |     +---- data
    |           |     +---- reference
    |           |     +---- run
    |           |     +---- src
    |           +---- dyn_vi_rhow_solver
    |                 +---- data
    |                 +---- reference
    |                 +---- run
    |                 +---- src
    +---- sysdep




# Quickest guide for the impatient #

The impatient may follow the lines below.

    $ cd $TOPDIR             # Replace with your extracted directory, such as ~/IcoAtmosBenchmark/
    $ export IAB_SYS=Linux64-intel-impi # Select from sysdep/
    $ cd sysdep
    $ edit Makedef.$IAB_SYS  # Check and edit for your environment.
    $ cd ../kernels/NICAM/
    $ make data              # for once only
    $ make                   # build executables
    $ make run               # run all executables


# Step by step procedure #

## Preparing input/validation data files ##

To run each kernel program and validate the results, you have to prepare
data file by downloading from our official site.

There is a shell script to download necessary files in `data/` under each
kernel program's directory, called `donwload.sh`.

*Note: Some of data files are more than 100MB large.*

Kernel communication does not need to download any files, and there are several .conf files.
You can change probrem size ('rlevel') by renaming one of `communication_gl04*.cnf` to `communication.cnf`.
Default is `gl04rl00`.

## Compile ##


You have to prepare `$(TOPDIR)/sysdep/Makedef.$(IAB_SYS)}` at first,
where `$(IAB_SYS)` is replaced by the system/complier
token. You have to define this as an environmental variable in your
shell.

There are several pre-defined files in `$(TOPDIR)/sysdep/` directory.
You can use one of them asis, or edit it to suit your environment.

Then you can compile each kernel program by using `make` command in each
`src/` directory.

*Note: This include file is used in common by all three models.*

## Run ##

There are several sample execution scripts in `run/` under each kernel
program's directory, such as `run-linux.sh`.

## Validation ##

Each kernel program outputs several lines to the standard output.
For example, the case of `dyn_diffusion` kernel are shown below.


    [KERNEL] dyn_diffusion
     *** Start  initialize
     *** Finish initialize
     *** Start  kernel
     ### Input ###
     +check[check_dscl      ] max=  6.1941315670898286E-08,min= -7.0374144752795210E-08,sum= -2.7407850230958588E-07
     +check[check_dscl_pl   ] max=  2.2139244811324830E-08,min= -6.0170678327656930E-11,sum=  1.8274579608905828E-07
     +check[scl             ] max=  1.6578626530298903E-11,min= -1.2670212860993856E-11,sum= -2.0289014286353776E-10
     +check[scl_pl          ] max=  1.9358849576453664E-11,min= -8.2853331106472899E-12,sum=  1.1445754720677440E-09
     +check[kh              ] max=  2.8341305529772246E+12,min=  6.7659597088284981E+10,sum=  6.0439053980501018E+17
     +check[kh_pl           ] max=  2.8334094314435532E+12,min=  6.7659597088284981E+10,sum=  4.3486317454839525E+14
     ### Output ###
     +check[dscl            ] max=  6.1941315670898286E-08,min= -7.0374144752795210E-08,sum= -2.7407850230958588E-07
     +check[dscl_pl         ] max=  2.2139244811324830E-08,min= -6.0170678327656930E-11,sum=  1.8274579608905828E-07
     ### Validation : point-by-point diff ###
     +check[check_dscl      ] max=  0.0000000000000000E+00,min=  0.0000000000000000E+00,sum=  0.0000000000000000E+00
     +check[check_dscl_pl   ] max=  0.0000000000000000E+00,min=  0.0000000000000000E+00,sum=  0.0000000000000000E+00
     *** Finish kernel
    
     *** Computational Time Report
     *** ID=001 : MAIN_dyn_diffusion               T=     0.028 N=      1
     *** ID=002 : OPRT_diffusion                   T=     0.028 N=      1


Check the lines below `### Validation : point-by-point diff ###` line,
those show difference between calculated output array and pre-calculated reference array.
These should be zero or enough small to be acceptable.

There are sample output log files in `reference/` in each kernel program
directory, for reference purpose.

