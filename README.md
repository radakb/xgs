XGS
===

eXpanded ensemble with Gibbs Sampling - A NAMD extension for expanded ensemble simulations

Description
===========

XGS implements expanded ensemble simulations for the popular software package NAnoscale Molecular Dynamics ([NAMD](http://www.ks.uiuc.edu/Research/namd/)). It utilizes NAMD's powerful Tcl interface to create new commands for these simulation types, generally with little added overhead. The package also includes calibration routines utilizing simple linear response relations which, in many cases, provide an excellent means for parameterizing a given expanded ensemble model.

Installation
============

XGS has two facets:
1) It is a generic library for expanded ensemble and Gibbs sampling routines.
2) It provides implementations of basic expanded ensemble simulation methods.

Neither of these require installation as such. In order to make use of provided methods, you need only "source" the relevant Tcl interface. These are meant to hide all of the "under the hood" functionality of XGS and only expose those parts that are directly needed for simulations. If you wish to implement your own bespoke expanded ensemble method, then you will need to source all or some of the XGS core namespace. Because NAMD essentially operates with its own Tcl interpreter, this must be done directly. One cannot, in general, make use of Tcl package management (e.g., "package require xxxx").

Examples
========
Coming Soon

Tests
=====
Coming Soon

Authors and Contributors
========================

* Brian Radak | brian.radak@gmail.com | brian.radak@anl.gov
