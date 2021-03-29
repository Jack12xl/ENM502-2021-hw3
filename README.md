#### Ling xie
#### ENM 502-001
#### 2021-03-08
## Assignment 3 Newton method and Arc-length continuation

### Introduction

In this assignment, we would solve a non-linear boundary-value problem defined on the unit-square domain 
$$
\begin{aligned}
&\begin{array}{l}
D=(0 \leq x \leq 1) \cup(0 \leq y \leq 1) \\
\nabla^{2} u+\lambda u(1+u)=0
\end{array}\\
&u(x, y)=0 \text { on all boundaries, }
\end{aligned}
$$


with **newton method**. We track the $\lambda$ from 
$$
0 \leq \lambda \leq 60
$$
with **analytic continuation(AyC)** and **arc-length continuation(ARCLC)**.



### Problem Setup and Formulation

##### Discretization and centered finite difference applied to this problem.

##### Problem Setup and Formulation

##### Methods to generate first two non-trivial solutions (analytical continuation)

##### Arc-Length Continuation

### This criterion is linked to a Learning Outcome Results and Discussion

From my circumstances, I notice that with different mixture of parameters( ARCLC step size, initial point, etc), the results vary rapidly.



##### Rush to bottom

| $\lambda$ iteration                    | L2norm of U $lambda$                          | U                                   | U                                   |
| -------------------------------------- | --------------------------------------------- | ----------------------------------- | ----------------------------------- |
| ![](./results/rush2bottom/lmbd_it.png) | ![](./results/rush2bottom/L2normU_lambda.png) | ![](./results/rush2bottom/fig1.png) | ![](./results/rush2bottom/fig2.png) |

Rush to top

| $\lambda$ iteration                 | L2norm of U $lambda$                       | U                                | U                                |
| ----------------------------------- | ------------------------------------------ | -------------------------------- | -------------------------------- |
| ![](./results/rush2top/lmbd_it.png) | ![](./results/rush2top/L2normU_lambda.png) | ![](./results/rush2top/fig1.png) | ![](./results/rush2top/fig2.png) |

Rush to top

| $\lambda$ iteration                  | L2norm of U $lambda$                        | U                                 | U                                 |
| ------------------------------------ | ------------------------------------------- | --------------------------------- | --------------------------------- |
| ![](./results/fluctuate/lmbd_it.png) | ![](./results/fluctuate/L2normU_lambda.png) | ![](./results/fluctuate/fig1.png) | ![](./results/fluctuate/fig2.png) |

Up and down

| $\lambda$ iteration               | L2norm of U $lambda$                     | U                              | U                              |
| --------------------------------- | ---------------------------------------- | ------------------------------ | ------------------------------ |
| ![](./results/updown/lmbd_it.png) | ![](./results/updown/L2normU_lambda.png) | ![](./results/updown/fig1.png) | ![](./results/updown/fig2.png) |



### Code

Code repository is [here](https://github.com/Jack12xl/ENM502-2021-hw3).