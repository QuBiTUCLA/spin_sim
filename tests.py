import unittest
import numpy as np

# class qgan_test(unittest.TestCase):
#     def test_build_qgan_v1(self):
#         """
#         Test that it GGAN object can return a random number.
#         """
#         train_set = weibull_dist(2., 2., 500, [0., 3.])
#         bounds = [0., 3.]
#         num_qubits = [2]
#         circuit_depth = 2
#         batch_size = 50
#         num_epochs = 20
#         myopt = ADAM(maxiter=10000, tol=1e-06, lr=0.001, beta_1=0.9, beta_2=0.99, noise_factor=1e-08, eps=1e-10,
#                      amsgrad=True)
#         my_qgan1 = build_qgan_v1(train_set, bounds, num_qubits, circuit_depth, batch_size, num_epochs, './')
#
#         self.assertIs(np.type(my_qgan1.random()), float, 'QGAN is working!')
#
#     def test_build_qgan_v2(self):
#         """
#         Test that it GGAN object can return a random number.
#         """
#         train_set = weibull_dist(2., 2., 500, [0., 3.])
#         bounds = [0., 3.]
#         num_qubits = [2]
#         circuit_depth = 2
#         batch_size = 50
#         num_epochs = 20
#         myopt = ADAM(maxiter=10000, tol=1e-06, lr=0.001, beta_1=0.9, beta_2=0.99, noise_factor=1e-08, eps=1e-10,
#                      amsgrad=True)
#         my_qgan2 = build_qgan_v2(train_set, bounds, myopt, num_qubits, circuit_depth, batch_size, num_epochs, './')
#
#         self.assertIs(np.type(my_qgan2.random()), float, 'QGAN is working!')


if __name__ == '__main__':
    unittest.main()