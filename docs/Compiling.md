


# Compile DINEOF

## Dynamic compilation

See `README.md`.

## Static compilation

Unfortunately, Ubuntu does no longer include static NetCDF libraries (https://bugs.launchpad.net/ubuntu/+source/netcdf/+bug/1698368)

```bash
sudo apt-get update
sudo apt-get install gfortran make libarpack2-dev libhdf5-dev
```

### NetCDF

```bash
wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-4.7.4.tar.gz
tar xvf netcdf-c-4.7.4.tar.gz 
cd netcdf-c-4.7.4/
./configure --disable-shared --enable-static  --disable-dap --prefix ~/opt/netcdf-4.7.4 LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib" CPPFLAGS="-I/usr/lib/x86_64-linux-gnu/hdf5/serial/include"
make
```

### NetCDF-Fortran

```bash
wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-fortran-4.5.2.tar.gz
tar -zxf netcdf-fortran-4.5.2.tar.gz
cd netcdf-fortran-4.5.2/
./configure --disable-shared --enable-static LDFLAGS="-L$HOME/opt/netcdf-4.7.4/lib/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib" CPPFLAGS="-I$HOME/opt/netcdf-4.7.4/include"  LIBS="-lhdf5_hl -lhdf5 -lm -lz" --prefix=$HOME/opt/netcdf-4.7.4
make test
make check
make install
```

### DINEOF

```bash
export PATH="$HOME/opt/netcdf-4.7.4/bin/:$PATH"
which nf-config
make STATIC=on
gfortran -static  -o dineof ReadMatrix.o ufileformat.o initfile.o stat.o norm.o dineof_utils.o smeanToZero.o smeanByRow.o svariExp.o ssvd_lancz.o dineof.o -larpack -llapack -lblas  -L$HOME/opt/netcdf-4.7.4/lib -lnetcdff -L$HOME/opt/netcdf-4.7.4/lib/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lnetcdf -lnetcdf -ldl -lm -lhdf5_hl -lhdf5 -lm -lz -lrt -pthread -ldl -lsz -laec
```

There is still a warning:
```/usr/bin/ld: /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.a(H5PLint.o): in function `H5PL__open':
(.text+0x4bd): Warnung: Using 'dlopen' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking```
