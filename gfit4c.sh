#!/bin/sh
# Compiles the fortran program 
f77 gfit4c.f -o gfit4c
# Executes the test run
gfit4c < gfit4c.inp > output
