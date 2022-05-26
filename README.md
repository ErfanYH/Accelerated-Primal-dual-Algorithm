# Accelerated-Primal-dual-Algorithm

This is a numerical example to showcase the performance of the Accelerated Primal-Dual algorithm with Backtracking (APDB) proposed in [https://arxiv.org/pdf/1803.01401.pdf](https://doi.org/10.1137/18M1213488). In this example, APDB algorithm is implemented for solving a Quadratic Constrained Quadratic Programming (QCQP). See section 5.1 of the published paper (or section 5.2 of the arXiv paper) for more details. 

APDB is proposed to solve general non-bilinear convex-concave saddle point problems. The proposed method achieves an optimal convergence rate of O(1/k) and O(1/k^2) in convex-concave and strongly convex-concave settings, respectively.


- APD_c corresponds to constant stepsizes (Theorem 2.1 part (I)).

- When the problem is strongly convex (mu>0) APD_sc uses non-constant stepsizes (Theorem 2.1 part (II)). 
