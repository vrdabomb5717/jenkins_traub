Title: Finding Roots of a Polynomial Using Jenkins-Traub  
Author: Varun Ravishankar  
Email: vr2263@columbia.edu  
Date: May 4th, 2013  

# Finding Roots of a Polynomial Using Jenkins-Traub
## by Varun Ravishankar

This program finds the roots of a polynomial using the Jenkins-Traub algorithm.


# Requirements

* Python 2.7
* Numpy

Both of these are available on the Columbia University CLIC lab machines.

# Usage

`$ python jenkins_traub.py -c [coefficients] -e [epsilon] -m [max_iterations]`

where

* `c` is a white-space delimited list of the coefficients. Coefficients can be complex.
* `e` is `EPSILON`, the epsilon value used to determine algorithm stopping conditions.
* `m` is `MAX_ITERATIONS`, the max number of iterations for the algorithm to run.

The coefficients can be complex. To specify a complex coefficient, use j instead of i. For example, to specify the number *1-i*, use `1-1j`. To specify *i*, use `1j`. If the complex part of a coefficient is 1, you must specify it. You may not use spaces to specify complex coefficients; `1 - 1j` will not work.

The max iterations by default is set to 1000, and epsilon by default is set to 1.0e-10.


# Examples

```
$ python jenkins_traub.py -c 1 -9.01 27.08 -41.19 32.22 -10.1
Done with Stage 1 after 5 iterations.
Done with Stage 2 after 2 iterations.
Done with Stage 3 after 7 iterations.
Done with Stage 1 after 5 iterations.
Done with Stage 2 after 2 iterations.
Done with Stage 3 after 6 iterations.
Done with Stage 1 after 5 iterations.
Done with Stage 2 after 2 iterations.
Done with Stage 3 after 6 iterations.
Done with Stage 1 after 5 iterations.
Done with Stage 2 after 2 iterations.
Done with Stage 3 after 4 iterations.
Roots: [(1.0000000020052087-4.6143098772373481e-10j), (1.0099999979899685+4.6254009123783083e-10j), (1.0000000000000868+0.99999999999002043j), (1.0000000000044402-0.99999999999106171j), (5.0000000000002958-6.7834626804597065e-14j)]

$ python -O jenkins_traub.py -c 1 -4 4
Done with Stage 1 after 5 iterations.
Done with Stage 2 after 2 iterations.
Done with Stage 3 after 1 iterations.
Roots: [(1.9999999999999996+5.5511151231257827e-16j), (2.0000000000000004-5.5511151231257827e-16j)]

$ python -O jenkins_traub.py -c 1+1j 10 0 2-5j 8 545
Done with Stage 1 after 5 iterations.
Done with Stage 2 after 4 iterations.
Done with Stage 3 after 3 iterations.
Done with Stage 1 after 5 iterations.
Done with Stage 2 after 4 iterations.
Done with Stage 3 after 3 iterations.
Done with Stage 1 after 5 iterations.
Done with Stage 2 after 4 iterations.
Done with Stage 3 after 4 iterations.
Done with Stage 1 after 5 iterations.
Done with Stage 2 after 8 iterations.
Done with Stage 3 after 4 iterations.
Roots: [(-2.1412696593208285+2.1444016657547995j), (-1.7159723186584845-2.0530205561017132j), (1.7466808485599237-1.7876397485309532j), (2.0691454602009416+1.763758058933768j), (-4.9585843307815534+4.9325005799440991j)]

```
