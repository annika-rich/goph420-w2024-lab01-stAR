# Lab Assignment 01 - Numerical Integration 
## GOPH420 Inversion and Parameter Estimation for Geophysicists
*Author:* A. Richardson

*Instructor:* B. Karchewski

*Semester:* Winter 2024

The purpose of this lab is to examine algorithm that perform numerical integration of functions for both discrete data and when the function of interest is known.

To use the code developed in this repository, refer to the requirements.txt file and install the listed package dependencies with the command:

    python -m pip install -r requirements.txt

The local goph420_lab01 package  will also need to be installed in developer mode in your local virtual environment. To install goph420_lab01 use the following command:

    pip install -e ./

The driver script for the functions developed in this repository can be found in examples/ directory. The functions script, integration.py, is located in the src/goph420_lab01/ directory. Figures generated by the different driver files can be found in figures/ directory. All unit tests and any figures generated by the tests are in the tests/ directory. Text files containing data used for any numerical integration algorithms can be found in the data/ directory.
