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


with **newton method**. We track the and iterate $\lambda$ from 
$$
0 \leq \lambda \leq 60
$$
with **analytic continuation(AyC)** and **arc-length continuation(ARCLC)**.



### Problem Setup and Formulation

##### Discretization and centered finite difference applied to this problem.

We discretize the `U` with finite difference method, in specifically a `30 X 30` uniform grid. 

##### Newton's Method

As a fixed-point iteration method, given a initial approximation at `U_0`, newton's method would approximate the solution in the derivative direction. 

The iteration goes like 
$$
\mathbf{u}^{k+1}=\mathbf{u}^{k}+\delta \mathbf{u}^{k}
$$
where $ \delta \mathbf{u}^{k} $ is calculated by 
$$
\left.{\mathbf{J}}\right|_{k}\left(\begin{array}{1}
{\partial \mathbf{u}}
\end{array}\right)_{k}=-\left({\partial {\mathbf{R}}}\right)_{k}
$$
where `J` is the Jacobin matrix.

For the initial value, when L2norm of $\lambda$ is near to zero, it could be visualized as a eigen value problem. The initial value near is like 



| lambda = 2 * pi^2            | lambda = 5 * pi^2            | lambda = 5 * pi^2              |
| ---------------------------- | ---------------------------- | ------------------------------ |
| ![](./results/init_2pi2.png) | ![](./results/init_5pi2.png) | ![](./results/init_5pi2_2.png) |



##### Methods to generate first two non-trivial solutions (analytical continuation)

We would use newton's method to solve for `U_0` from first initial guess(see image above). To calculate `U_1` , we use **AyC** to calculate the initial guess `U_1_guess`, then use newton's method to solve the solution at `lambda_1`.

##### Arc-Length Continuation

During iterating the `lambda`, the `J ` would become nearly singular during turning point. To prevent this, we introduce a new independent variable `s`.
$$
(\delta S)^{2}=(\delta \lambda)^{2}+\|\delta \mathbf{u}\|_{2}^{2}
$$
and start to iterate on `s` for the rest steps.

### This criterion is linked to a Learning Outcome Results and Discussion

From my circumstances, I notice that with different mixture of parameters( **ARCLC** step size, initial point, etc), the results vary rapidly.

#### Failed cases

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

#### Relatively make sense cases

| $\lambda$ iteration                 | L2norm of U $lambda$                       | U                                | U                                |
| ----------------------------------- | ------------------------------------------ | -------------------------------- | -------------------------------- |
| ![](./results/can_flip/lmbd_it.png) | ![](./results/can_flip/L2normU_lambda.png) | ![](./results/can_flip/fig1.png) | ![](./results/can_flip/fig2.png) |

| $\lambda$ iteration              | L2norm of U $lambda$                    | U                             | U                             |
| -------------------------------- | --------------------------------------- | ----------------------------- | ----------------------------- |
| ![](./results/2pi_2/lmbd_it.png) | ![](./results/2pi_2/L2normU_lambda.png) | ![](./results/2pi_2/fig1.png) | ![](./results/2pi_2/fig2.png) |

### Conclusion

This is the expected whole graph of `L2 norm of U` vs `iteration`.

![expected_results](./results/expected_results.jpg)

In my current version or ARCLC, the full newton's method(with `J_hat`) could hardly converge. I think that's the reason why the `lambda` vs iteration is not stable and smooth.

On the other hand, we manage to see that the `U` could manage to switch on another branch(see the last two examples.).  

For example, in the last example, the results jumps from the red branch to the green branch.

In some circumstances(the first results), the `U` could sometimes move to `2pi^2` branch. 

### Code

Code repository is [here](https://github.com/Jack12xl/ENM502-2021-hw3).