# What is this ? #

**IcoAtmosBenchmark** is the benchmark suite of compulational kernels from
three icosahedral atmospheric model developed and used in world wide,
**NICAM**, **DYNAMICO** and **ICON**.

Kernels are extracted to be individual executables and you are able to
run them separately.

# kernels #

In this package, kernels from **DYNAMICO** are as belows.
See the description manual for the details of each kernel.

* `comp_pvort`
* `comp_geopot`
* `comp_caldyn_horiz`
* `comp_caldyn_vert`

# Directory tree #

Kernel programs are stored in individual directories,
under each model's directory tree, such as
 `${TOPDIR}/kernels/DYNAMICO/comp_pvort`.

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
    |     |     +---- comp_pvort
    |     |     |     +---- conf
    |     |     |     +---- reference
    |     |     |     +---- run
    |     |     |     +---- src
    |     |     +---- comp_geopot
    |     |     |     +---- conf
    |     |     |     +---- reference
    |     |     |     +---- run
    |     |     |     +---- src
    |     |     +---- comp_caldyn_horiz
    |     |     |     +---- conf
    |     |     |     +---- reference
    |     |     |     +---- run
    |     |     |     +---- src
    |     |     +---- comp_caldyn_vert
    |     |           +---- conf
    |     |           +---- reference
    |     |           +---- run
    |     |           +---- src
    |     |
    |     +---- ICON
    |     |     +---- .....
    |     |     +---- .....
    |     +---- NICAM
    |           +---- .....
    |           +---- .....
    +---- sysdep



# Quickest guide for the impatient #

The impatient may follow the lines below.

    $ cd $TOPDIR             # Replace with your extracted directory, such as ~/IcoAtmosBenchmark/
    $ export IAB_SYS=Linux64-intel-impi # Select from sysdep/
    $ cd sysdep
    $ edit Makedef.$IAB_SYS  # Check and edit for your environment.
    $ cd ../kernels/DYNAMICO
    $ make data              # for once only
    $ make                   # build executables
    $ make run               # run all executables


# Step by step procedure #

## Preparing input/validation data files ##

To run each kernel program and validate the results, you have to prepare
data file by downloading from our official site.

There is a shell script to download necessary files in `data/` under each
kernel program's directory, called `donwload.sh`.

*Note: Some of data files are several MB large.*


## Compile ##


You have to prepare `$(TOPDIR)/sysdep/Makedef.$(IAB_SYS)}` at first,
where `$(IAB_SYS)` is replaced by the system/complier
token. You have to define this as an environmental variable in your
shell.

There are several pre-defined files in `$(TOPDIR)/sysdep/` directory.
You can use one of them asis, or edit it to suit your environment.

Then you can compile each kernel program by using `make` command in each
`src/` directory.

*Noth: This include file is used in common by all three models.*

## Run ##

There are several sample execution scripts in `run/` under each kernel
program's directory, such as `run-linux.sh`.

## Validation ##

Each kernel program outputs several lines to the standard output.
For example, the case of `comp_pvort` kernel are shown below.


    [KERNEL] comp_pvort
     *** Start  initialize
                    iim, jjm, llm:    23    25    19
                 ij_begin, ij_end:    48   528
         ij_begin_ext, ij_end_ext:    24   552
                 ll_begin, ll_end:     1    19
            t_right, t_rup, t_lup:     1    23    22
         t_left, t_ldown, t_rdown:    -1   -23   -22
            u_right, u_rup, u_lup:     0  1173   575
         u_left, u_ldown, u_rdown:    -1  1150   553
               z_rup, z_up, z_lup:   598     0   597
         z_ldown, z_down, z_rdown:   -23   575   -22
                       caldyn_eta:     1
                                g:     9.80000000
     +check[Av              ] max=  4.1228713627140027E+11,min=  0.0000000000000000E+00,sum=  3.2428753277257527E+13
     +check[de              ] max=  4.5171816240714993E+06,min=  0.0000000000000000E+00,sum=  4.7785815753077912E+08
     +check[Riv2            ] max=  3.8193271158709069E-01,min=  0.0000000000000000E+00,sum=  9.6499999999999977E+02
     +check[fv              ] max=  8.2383275804860789E-05,min= -6.6879410680009186E-05,sum=  4.1115175112148659E-03
     +check[mass_dak        ] max=  3.9864758943842335E+03,min= -4.2200064724847380E+03,sum=  1.1368683772161603E-13
     +check[mass_dbk        ] max=  1.6280745944237918E-01,min=  0.0000000000000000E+00,sum=  1.0000000000000000E+00
     *** Finish initialize
     *** Start kernel
     ### check point iteration:        1000
     ### Input ###
     +check[u               ] max=  1.0215854276981204E+01,min= -2.0878929653719833E+00,sum=  5.0175693885050535E+01
     +check[ps              ] max=  1.0000000000000000E+05,min=  1.0000000000000000E+05,sum=  5.7500000000000000E+07
     +check[theta_rhodz     ] max=  3.9393370687019045E+05,min=  0.0000000000000000E+00,sum=  1.8099621340626464E+09
     +check[theta_prev      ] max=  8.0139914420291746E+02,min=  0.0000000000000000E+00,sum=  3.8582633571973117E+06
     +check[rhodz_prev      ] max=  1.2306877011993038E+03,min=  0.0000000000000000E+00,sum=  5.3979591836733194E+06
     +check[qu_prev         ] max=  1.0339537867296609E-06,min= -8.4408169682701225E-07,sum=  3.9419811615778674E-04
     +check[qv_prev         ] max=  1.0397552030841796E-06,min= -8.4408169685381862E-07,sum=  2.6984926372303133E-04
     ### Output ###
     +check[theta           ] max=  8.0139914420291746E+02,min=  0.0000000000000000E+00,sum=  3.8582633571973117E+06
     +check[rhodz           ] max=  1.2306877011993038E+03,min=  0.0000000000000000E+00,sum=  5.3979591836733194E+06
     +check[qu              ] max=  1.0626772908333491E-06,min= -8.5290650439776975E-07,sum=  3.9932412024562623E-04
     +check[qv              ] max=  1.0864352277577492E-06,min= -8.9078811385791910E-07,sum=  2.7528879622919179E-04
     ### final iteration:        1000
     ### Validation : grid-by-grid diff ###
     +check[theta           ] max=  0.0000000000000000E+00,min=  0.0000000000000000E+00,sum=  0.0000000000000000E+00
     +check[rhodz           ] max=  0.0000000000000000E+00,min=  0.0000000000000000E+00,sum=  0.0000000000000000E+00
     +check[qu              ] max=  0.0000000000000000E+00,min=  0.0000000000000000E+00,sum=  0.0000000000000000E+00
     +check[qv              ] max=  0.0000000000000000E+00,min=  0.0000000000000000E+00,sum=  0.0000000000000000E+00
     *** Finish kernel
    
     *** Computational Time Report
     *** ID=001 : MAIN_comp_pvort                  T=     0.232 N=   1000


Check the lines below `### Validation : grid-by-grid diff ###` line,
those show difference between calculated output array and pre-calculated reference array.
These should be zero or enough small to be acceptable.

There are sample output log files in `reference/` in each kernel program
directory, for reference purpose.
