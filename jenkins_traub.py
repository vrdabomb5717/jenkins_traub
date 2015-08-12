#!/usr/bin/env python

"""Find the roots of a polynomial using the Jenkins-Traub algorithm."""


from __future__ import division
from __future__ import print_function

import argparse
import re
from functools import partial

import numpy as np


class ConvergenceError(Exception):
    pass


class RootFound(Exception):
    pass


def evaluate(coefficients):
    """Evaluate a polynomial with given coefficients in highest to lowest order."""
    return lambda s: synthetic_division(coefficients, s)[1]


def newton(x_0, function, first_derivative, epsilon, max_iterations):
    """Find the root of a polynomial using Newton-Raphson iteration."""
    x_i = x_0
    num_iterations = 0
    x_i_next = x_i - (function(x_i) / first_derivative(x_i))
    num_iterations += 1

    while abs(x_i_next - x_i) > epsilon and num_iterations < max_iterations:
        x_i = x_i_next
        x_i_next = x_i - (function(x_i) / first_derivative(x_i))
        num_iterations += 1

    return x_i


def synthetic_division(coefficients, s):
    """Perform synthetic division and evaluate a polynomial at s."""
    deflated = np.empty(len(coefficients) - 1, dtype=np.complex128)

    deflated[0] = coefficients[0]

    for i, a_i in enumerate(coefficients[1:-1], start=1):
        deflated[i] = a_i + deflated[i - 1] * s

    evaluation = coefficients[-1] + deflated[-1] * s

    return deflated, evaluation


def next_step(a, H_bar_lambda, s, epsilon, generate_t=True):
    """Generate the next H_bar_lambda and t, if desired."""
    p, p_at_s = synthetic_division(a, s)
    h_bar, h_bar_at_s = synthetic_division(H_bar_lambda, s)

    # If we found a root, short circuit the other logic.
    # Used for when we found a root exactly, and we'd end up dividing
    # by 0 when finding s_lambda
    if np.absolute(p_at_s) < epsilon:
        raise RootFound()

    if np.absolute(h_bar_at_s) < epsilon:
        h_bar_at_s += epsilon / 100

    t = None
    if generate_t:
        t = s - p_at_s / h_bar_at_s

    h_bar = np.insert(h_bar, 0, 0)  # watch for polynomial length differences
    return p - (p_at_s / h_bar_at_s) * h_bar, t


def stage1(a, H_lambda, epsilon, max_iterations):
    """Perform Stage 1, the no-shift process of Jenkins-Traub.

    Returns the deflated polynomial, and the number of iterations the stage took.

    """
    H_bar_lambda = H_lambda / H_lambda[0]

    for num_iterations in xrange(max_iterations):
        try:
            H_bar_lambda, _ = next_step(a, H_bar_lambda, 0, epsilon, False)
        except RootFound:
            break

    return H_bar_lambda, num_iterations


def stage2(a, H_lambda, s, epsilon, max_iterations):
    """Perform Stage 2, the fixed shift process of Jenkins-Traub.

    Returns the deflated polynomial, and the number of iterations the stage took.

    """
    t = t_prev = float('inf')
    num_iterations = 0

    H_bar_lambda = H_lambda / H_lambda[0]

    while True:
        t_prev, t_prev_prev = t, t_prev

        try:
            H_bar_lambda, t = next_step(a, H_bar_lambda, s, epsilon)
        except RootFound:
            break

        num_iterations += 1

        condition1 = np.absolute(t_prev - t_prev_prev) <= 0.5 * np.absolute(t_prev_prev)
        condition2 = np.absolute(t - t_prev) <= 0.5 * np.absolute(t_prev)
        condition3 = num_iterations > max_iterations

        if (condition1 and condition2) or condition3:
            break

    return H_bar_lambda, num_iterations


def stage3(a, H_L, s, epsilon, max_iterations):
    """Perform Stage 3, the variable-shift recurrence of Jenkins-Traub.

    Returns the deflated polynomial, the root, and the number of
    iterations the stage took.

    """
    polynomial = evaluate(a)
    H_bar_coefficients = H_L / H_L[0]
    H_bar = evaluate(H_bar_coefficients)
    s_L = s - polynomial(s) / H_bar(s)

    H_bar_lambda = H_bar_coefficients.copy()
    s_lambda = s_L
    s_lambda_prev = complex('inf')

    num_iterations = 0

    while np.absolute(s_lambda - s_lambda_prev) > epsilon and num_iterations < max_iterations:
        p, p_at_s_lambda = synthetic_division(a, s_lambda)
        h_bar, h_bar_at_s_lambda = synthetic_division(H_bar_lambda, s_lambda)
        h_bar = np.insert(h_bar, 0, 0)  # watch for polynomial length differences

        # If we found a root, short circuit the other logic.
        # Used for when we found a root exactly, and we'd end up dividing
        # by 0 when finding s_lambda
        if np.absolute(p_at_s_lambda) < epsilon:
            return H_bar_lambda, s_lambda, num_iterations

        H_bar_lambda_next = p - (p_at_s_lambda / h_bar_at_s_lambda) * h_bar
        H_bar_next = evaluate(H_bar_lambda_next)

        num_iterations += 1

        s_lambda, s_lambda_prev = (s_lambda - p_at_s_lambda / H_bar_next(s_lambda)), s_lambda
        H_bar_lambda = H_bar_lambda_next

    if num_iterations >= max_iterations:
        print('Stage 3 could not converge after {} iterations'.format(num_iterations))
        raise ConvergenceError(num_iterations)

    return H_bar_lambda, s_lambda, num_iterations


def jenkins_traub_inner(coefficients, epsilon, max_iterations, do_stage_one=True):
    """Find the smallest root of a polynomial using the Jenkins-Traub algorithm.

    Returns the deflated polynomial and the (approximately) smallest root.

    """
    # Find zero root, if any
    if len(coefficients) > 0 and coefficients[-1] == 0:
        s_lambda = 0
        return coefficients, s_lambda

    if len(coefficients) < 2:
        return

    a = np.array(coefficients)

    powers = np.array(xrange(len(coefficients) - 1, -1, -1))

    # Scale coefficients so leading coefficient is 1
    if a[0] != 1:
        a /= a[0]

    if len(a) == 2:
        s_lambda = -a[-1] / a[0]
        return None, s_lambda

    # H^0 is the derivative of the polynomial
    H_0 = powers * a
    H_0 = H_0[:-1]

    H_lambda = H_0.copy()
    s_lambda = None

    # At this point, you could perform Stage 1 to accentuate small zeros.
    # We choose to do this and then move on to Stage 2.
    if do_stage_one:
        H_lambda, num_iterations1 = stage1(a, H_0, epsilon, 5)
        print("Done with Stage 1 after {} iterations.".format(num_iterations1))

    # Get s using a modified Cauchy polynomial and Newton's iteration.
    # Modified polynomial has coefficients that are the moduli of the
    # original polynomial, but the last root is -1.
    modified_coefficients = np.absolute(a)
    modified_coefficients[-1] *= -1
    modified_derivative = powers * modified_coefficients
    modified_derivative = modified_derivative[:-1]

    mod_function = evaluate(modified_coefficients)
    mod_derivative = evaluate(modified_derivative)

    # Only get beta up to two decimal places.
    x_0 = 1
    beta = newton(x_0, mod_function, mod_derivative, 0.01, 500)

    while True:
        try:
            phi_random = complex(np.random.random() * 2 * np.pi)
            s_lambda = np.exp(phi_random * (0+1j)) * beta
            # s_lambda = 10000+10000j  # use for testing bad convergence
            H_lambda, num_iterations2 = stage2(a, H_lambda, s_lambda, epsilon, 100)
            print("Done with Stage 2 after {} iterations.".format(num_iterations2))

            H_lambda, s_lambda, num_iterations3 = stage3(a, H_lambda, s_lambda, epsilon, max_iterations)
            print("Done with Stage 3 after {} iterations.".format(num_iterations3))
            return H_lambda, s_lambda
        except ConvergenceError:
            max_iterations *= 2
            continue


def jenkins_traub(coefficients, epsilon, max_iterations):
    """Find the roots of a polynomial using the Jenkins-Traub algorithm.

    Wrapper method for the Jenkins-Traub algorithm. Loops through the
    coefficients, finds a root, divides the coefficients by
    that smallest root to find a new set of coefficients, and then starts over.

    Returns a list of the roots.

    """
    # Remove leading coefficients equal to 0
    while len(coefficients) > 0 and coefficients[0] == 0:
        coefficients.pop(0)

    if len(coefficients) == 0:
        print("Fed the zero polynomial. Roots are all complex numbers.")
        return
    elif len(coefficients) == 1:
        print("Fed a constant. There are no roots.")
        return

    for x in xrange(len(coefficients) - 1):
        if len(coefficients) > 1:
            _, s = jenkins_traub_inner(coefficients, epsilon, max_iterations)
            coefficients, _ = synthetic_division(coefficients, s)
            yield s


def parse_args():
    def positive(data, datatype=float):
        """Check if a string converts to a positive number."""
        msg = "{} is not a positive number.".format(data)
        try:
            value = datatype(data)
        except ValueError:
            raise argparse.ArgumentTypeError(msg)
        if value < 0:
            raise argparse.ArgumentTypeError(msg)
        return value

    positive_float = partial(positive, datatype=float)
    positive_int = partial(positive, datatype=int)

    parser = argparse.ArgumentParser(description="""Find the roots of a polynomial
                                     using the Jenkins-Traub algorithm.""")
    parser._negative_number_matcher = re.compile(r'^-.+$')
    parser.add_argument("-c", type=complex, dest='coefficients', nargs='+',
                        default=[1, -9.01, 27.08, -41.19, 32.22, -10.1],
                        # default=[1, -4, 4],
                        # default=[1, 0, 0, 0, -16],
                        help="The coefficients of the polynomial")
    parser.add_argument("-e", type=positive_float, dest='epsilon', default=1.0e-10,
                        help="Epsilon value used to determine algorithm stopping conditions.")
    parser.add_argument("-m", type=positive_int, dest='max_iterations', default=1000,
                        help="The max number of iterations for the algorithm to run.")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    coefficients = args.coefficients
    epsilon = args.epsilon
    max_iterations = args.max_iterations

    solution = jenkins_traub(coefficients, epsilon, max_iterations)
    roots = [x for x in solution]

    print("Roots: {}".format(roots))


if __name__ == '__main__':
    main()
