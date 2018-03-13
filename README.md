P2.6 - The Finite Element Method Using deal.II
===============================================
## Lecturers: Jean-Paul Pelteret and Luca Heltai

This repository contains the assignements and workspaces for the
course P2.6

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


