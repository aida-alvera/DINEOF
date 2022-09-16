## Troubleshooting

### "Segmentation fault"

When compiling with `ifort`, you may encounter a "segmentation fault" error when running DINEOF. Something like this:
```bash
****************************************************
 Numerical data read
 You entered the values:
 number of EOF modes you want to compute           2
 maximal size for the Krylov subspace           8

 You asked not to  normalise of the input matrices

 The right and left EOFs will be written in directory Output/

****************************************************
You entered filenames seacoos2005.avhrr
										 seacoos2005.avhrr.mask
Data error in UINQUIRE, not a conform file: "seacoos2005.avhrr".
Segmentation fault
```

#### Solution
The problem appears to be that `ifort` allocates too few stack memory by default. If you type the following command in a terminal:
```bash
ulimit -s 1000000
```
That may solve the problem (the number needs to be a very high number, whatever you choose)

### "Error = -3"

Another common problem is when you start to see this message appear in the screen when you run DINEOF:
```bash
Error with _saupd, info =           -3
Check documentation in _saupd
```bash

#### Solution
This one is also simple to solve: you need to specify in `dineof.init` the number of modes to be calculated (`nev`), and the _Kryvlov_ space (`ncv`), which need to be smaller than the temporal size of your data, AND, `ncv > nev+5`. If one or both of these conditions is not met, the above message will appear.

### "Error = -8"

A less common problem:
```bash
Error with _saupd, info =           -8
Check documentation in _saupd
```

This problem is related to the filtering of the time covariance matrix: if the time vector is not well written (i.e., the values are not all distinct), the filtering cannot be performed properly. As the time is in simple precision, even if the values are distinct, they can be read as identical values.

Here is an example of a time vector which lead to the problem:
```bash
ncdump -v time bstot_mumm.nc
data:
time = 731582, 731582.020833333, 731582.041666667, 731582.0625,
```

#### Solutions

* Deactivate the filtering (`B_DIFF` option in file `ppdef.h`)
* Write the time in double precision

### netCDF compilation error

In some systems, the `Fortran` part of netCDF is in a separate netCDF library, called "netcdff". In those, it is needed to add `-lnetcdff` (two f's) before the `-lnetcdf` option, in the compilation command:
```bash
-lnetcdff -lnetcdf
```
If you use the DINEOF Makefile, this means to add this option to the "DINEOF_LIBRARIES" variable:
```bash
DINEOF_LIBRARIES =  -L/usr/lib -larpack  -llapack  -lblas  -lnetcdff -lnetcdf
```
of the compilation file (`Linux-gfortran.mk`, etc)

### Compilation error due to ARPACK

If, when compiling `DINEOF`, you find the following error:
```bash
/home/user/ARPACK/libarpack.a(second.o): In function 'second_':
second.f:(.text+0x16): undefined reference to 'etime_'
collect2: ld returned 1 exit status
make: *** [dineof] Error 1
```
you will need to edit file `'second.f'` located in your `$HOME/ARPACK` folder. There, you should comment the following line by adding an asterisk at the beginning:
```bash
EXTERNAL           ETIME
```
becomes:
```bash
* EXTERNAL           ETIME
```

### Only missing values with netCDF files

When using netCDF files, please use the missing_value or `_FillValue` with a real number (float) and not with NaN (though this is now tolerated in netCDF files it is not yet implemented in DINEOF and can result in interpretations of your file having no data)


### Integer overflow when running the small example

You get the following error when you run the small example.
```bash
Fortran runtime error: Integer overflow when calculating the amount of memory to allocate.
```
Note that the data files of the small example are in "big endian". You have to leave the format to "big_endian" in the Makefile to run the small example. You can only use another format when you create yourself the data set.


###  Compilation fails with `could not open netcdf.inc`

Make sure to install the *NetCDF Fortran development package*. It is called `libnetcdff-dev` on Debian/Ubuntu (note the 2nd `f`) or `netcdf-fortran-devel` on Fedora/Red Hat/Centos/OpenSuse .

For example, you can install it with the following command on Debian/Ubuntu:

```bash
sudo apt install libnetcdff-dev
```
