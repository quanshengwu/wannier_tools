Author: QuanSheng Wu (wuquansheng@gmail.com)
License: GPL V3

This is a code for twisted graphene systems using $k\cdot p$ model. It can be used to generate results for this paper
"ShengNan Zhang, Bo Xie, QuanSheng Wu, Jianpeng Liu, Oleg V. Yazyev, Chiral Decomposition of Twisted Graphene Multilayers with Arbitrary Stacking, arXiv:2012.11964 (2020)"

1. Installation
This code is written in Fortran, Message Passing Interface (MPI) is supported. Before the compilation, 
You need to modify the Makefile file to set up the Fortran compiler, Lapack and Blas libraries. then type

```
make
```

After the compilation, an excutable file ***tg_kpgen*** will generated.

2. The control parameters are included in the file *system.in*. 

3. run the code

```
./tg-kpgen
```
