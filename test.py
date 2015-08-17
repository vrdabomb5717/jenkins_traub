"""Use numpy to verify the algorithm's results.

$ python2 -m unittest test
"""
from unittest import TestCase

import numpy as np

import jenkins_traub
# silence the print statements
setattr(jenkins_traub, 'print', lambda _: None)


class TestJenkinsTraub(TestCase):
    def jenkins_traub(self,
                      coefficients, epsilon=1.0e-10, max_iterations=1000):
        return jenkins_traub.jenkins_traub(
            coefficients, epsilon, max_iterations
        )

    def check_results(self, coefficients):
        """Fail if the roots are wrong.

        Ideally there would be something like assertItemsAlmostEqual, but
        there isn't.

            self.assertItemsAlmostEqual(jt, roots)

        Even worse, sorting the arrays and then comparing elementwise doesn't
        work because of the complex numbers.

            np.testing.assert_almost_equal(np.sort(jt), np.sort(roots)))
        """
        jt = list(self.jenkins_traub(coefficients))
        roots = np.roots(coefficients)

        self.assertEqual(len(jt), len(roots))

        # Iterate through the Jenkins-Traub roots, removing roots from the
        # numpy roots if there is a match. If no match is found, the test
        # fails.
        expected_roots = np.copy(roots)
        for jt_root in jt:
            for index, root in enumerate(expected_roots):
                if np.isclose(jt_root, root):
                    expected_roots = np.delete(expected_roots, index)
                    break
            else:
                self.fail('roots differ: {} != {}'.format(jt, roots))

    def test_case_one(self):
        c = [1, -9.01, 27.08, -41.19, 32.22, -10.1]
        self.check_results(c)

    def test_case_two(self):
        c = [1, -4, 4]
        self.check_results(c)

    def test_case_three(self):
        c = [1+1j, 10, 0, 2-5j, 8, 545]
        self.check_results(c)
