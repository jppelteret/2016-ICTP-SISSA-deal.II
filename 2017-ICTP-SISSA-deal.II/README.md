P2.13 - The Finite Element Method Using deal.II
===============================================

This repository contains the assignements and workspaces for the
course P2.13

New material will be uploaded frequently,

```
  git remote add P2.13_seed https://github.com/sissa/P2.13_seed.git
```

or


```
  git remote add P2.13_seed git@github.com:sissa/P2.13_seed.git
```

and to

```
  git pull P2.13_seed master
```

in the morning.

Running deal.II on Ulysses
==========================

If you have access to Uylsses, you can add the following to your `.bashrc`:

	. /home/mathlab/candi-gnu-6.2.0/configuration/enable.sh

This will export all libraries required by `deal.II`, and the latest stable  version of `deal.II` (v8.5.0).

Running deal.II using docker
============================

```
docker run -t -i dealii/dealii:v8.5.0-gcc-mpi-fulldepscandi-debugrelease

```
