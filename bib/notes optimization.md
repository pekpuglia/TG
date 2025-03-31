# Notes

## Questions

* How to show a problem is convex or not?
* How to show configurations where a method is better than the other for high dimensional configurations?
* Switching surface?
* primer vector?
* T_final < T_max?
* Pontryagin's Minimum Principle
* solve sequence of problems starting with impulse + Lambert and half the thrust until desired thrust is obtained

## Bruce Conway 1.2.1

$\dot x = f(x, u, t)$

$x_i(0)$ given

$\Psi[x(T), T] = 0$ terminal condition

Objective:

$J = \phi[x(T), T] + \int_o^T L[x, u, t] dt$

### Necessary conditions for extrema

(Bryson, Ho)
Hamiltonian:

$H = L + \lambda^Tf$

Necessary conditions:
$\begin{cases} \dot \lambda = - \frac{\partial H}{\partial x} \\ \lambda(T) = [\frac{\partial \phi}{\partial x} + \nu^T \frac{\partial \Psi}{\partial x}]_{t=T} \\ \frac{\partial H}{\partial u} = 0 \end{cases}$

where

* $\lambda(t)$ is a costate trajectory
* $\nu$ is a multiplier vector TBD

**OBS**: choose $L = \Gamma$ (acceleration) $\rightarrow$ cost = $\Delta V$. In this case the Hamiltonian is linear in u so $\frac{\partial H}{\partial u} = 0$ is never true. Use Pontryagin's Minimum Principle. Therefore $\hat u \parallel - \lambda_v$ (PRIMER VECTOR). 

$p(t) = - \lambda_v(t)$

Primer vector equations:

$\begin{cases} \ddot p = \frac{\partial g(r)}{\partial r} p = G(r) p \\ p(t_f) = -\nu^T \frac{\partial \Psi}{\partial v(t_f)} \\ \dot p(t_f) = \nu^T \frac{\partial \Psi}{\partial r(t_f)} \end{cases}$

$\hat u = p(t) / \lVert p(t) \rVert$

This leads to switching function:

$S(t) = \lVert p \rVert - 1$

Bang-bang:

$\Gamma = \begin{cases} \Gamma_{max}, S > 0 \\ 0, S < 0 \end{cases}$

On impulsive trajectory, $\lVert p \rVert \leq 1$.

* Improve non opt traj with primer vector => Lion and Handelsman

* H not explicitly time depedent => should be constant (??), verify in numerical solution

* shooting, finite difference, collocation
* problem: determine initial costate???
* static dynamic control: Mystic

### Discretization

common, high dimensional NLP

### Evolutionary Algorithms

no initial guess, discrete

## Bruce Conway Primer Vector

Minimum propellant consumption:

$\max m(T)$

### Constant specific impulse 

Through $\Delta V(m_f)$:

$\min J = \min \int_0^T \Gamma(t) dt$



### Variable Specific Impulse

solar panels, radioisotope, ...

Power limited (<https://www.sciencedirect.com/science/article/pii/S009457651200166X>, <https://www.sciencedirect.com/science/article/pii/0094576589901069> **ITA**)

$\arg \max m(t_f) = \arg \min \int_0^{t_f} \Gamma^2 dt$

Optimal: $\Gamma(t) = p(t)$



**Question** jacobi no-conjugate-point condition page 23?

### Primer vector equation

$\ddot p = \frac{\partial g(r)}{\partial r} p = G(r) p$

Transition matrix $p(t) = \Phi(t, t_0) p(t_0)$ for costate and state in coasting arcs.

**Question** GLANDORF????

$\Phi^{-1}(t, t_0) = - J \Phi^T(t, t_0) J$

$J = \begin{bmatrix} 0_3 & I_3 \\ -I_3 & 0_3 \end{bmatrix}$

### Application of Primer Vector Theory to an Optimal Impulsive Trajectory

terminal coast: optimal transfer time < specified transfer time

* check sign $\dot p$ at extrema to see if initial/final coast are needed.
* if p > 1 anywhere, an extra impulse is needed