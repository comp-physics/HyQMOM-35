# Mathematical Foundations

This section provides the mathematical background for the Hyperbolic Quadrature Method of Moments (HyQMOM) as implemented in HyQMOM.jl.

## Kinetic Equation and Moment Systems

### 3D Kinetic Equation

HyQMOM.jl solves the 3D kinetic equation for the velocity distribution function (VDF) $f(t,\mathbf{x},\mathbf{u})$:

```math
\partial_t f + \mathbf{u} \cdot \nabla_\mathbf{x} f = \frac{1}{\text{Kn}} Q(f)
```

where:
- $\mathbf{u} = (u_1, u_2, u_3)$ is the velocity vector
- $\mathbf{x} = (x_1, x_2, x_3)$ is the spatial position
- $\text{Kn}$ is the Knudsen number
- $Q(f)$ is the collision operator (BGK model)

### Moment Definition

The velocity moments are defined as:
```math
m_{ijk} = \int_{\mathbb{R}^3} u_1^i u_2^j u_3^k f \, du_1 du_2 du_3
```

For the 35-moment system, we use all moments up to fourth order: $i+j+k \leq 4$.

### Central and Standardized Moments

**Central moments** remove the mean velocity:
```math
c_{ijk} = \int_{\mathbb{R}^3} (u_1-\bar{u}_1)^i (u_2-\bar{u}_2)^j (u_3-\bar{u}_3)^k f \, du_1 du_2 du_3
```

**Standardized moments** normalize by the variance:
```math
s_{ijk} = \frac{c_{ijk}}{c_{200}^{i/2} c_{020}^{j/2} c_{002}^{k/2}}
```

## HyQMOM Closure Methodology

### Orthogonal Polynomials

HyQMOM uses orthogonal polynomials $Q_n(x)$ constructed from the standardized moments via the three-term recurrence:

```math
Q_{n+1}(x) = (x - a_n) Q_n(x) - b_n Q_{n-1}(x)
```

The recurrence coefficients are:
```math
a_n = \frac{\langle x Q_n^2 \rangle}{\langle Q_n^2 \rangle}, \quad b_n = \frac{\langle Q_n^2 \rangle}{\langle Q_{n-1}^2 \rangle}
```

### Multivariate Closure

For multidimensional systems, closures are found using orthogonality:
```math
\langle \mathcal{U}_i \mathcal{V}_j \mathcal{W}_k \rangle = 0 \quad \text{for } i+j+k = 2n+1
```

where $\mathcal{U}_i$, $\mathcal{V}_j$, $\mathcal{W}_k$ are orthogonal polynomials in each direction.

### Key Properties

**Global Hyperbolicity**: All eigenvalues of the flux Jacobian are real

**Moment Realizability**: Solution moments correspond to a valid probability distribution

**Exact for Separable VDFs**: Closures are exact when $f = f_1(u_1)f_2(u_2)f_3(u_3)$

## Moment Realizability

### Realizability Constraints

For moments to be realizable (correspond to a valid VDF), certain matrix determinants must be positive:

**2nd-order constraint:**
```math
|\Delta_1| = 1 + 2s_{110}s_{101}s_{011} - (s_{110}^2 + s_{101}^2 + s_{011}^2) > 0
```

**4th-order constraint:**
```math
|\Delta_2| = |\Delta_1| \cdot |\Delta_2^*| > 0
```

where $\Delta_2^*$ is a $6 \times 6$ matrix involving all moments up to 4th order.

### Correction Algorithm

When non-realizable moments are detected:

**Check 2nd-order moments**: Ensure $|\Delta_1| > 0$

**Classify boundary type**: Corner, edge, or face of moment space

**Apply corrections**: Scale or reset moments to boundary values

**Verify 4th-order realizability**: Check $|\Delta_2^*| > 0$

### Physical Interpretation

**Interior of moment space**: General non-equilibrium distributions

**Boundary of moment space**: Highly non-equilibrium states (e.g., crossing jets)

**Corners**: Perfectly correlated velocities, $s_{110} = s_{101} = s_{011} = \pm 1$

## Hyperbolicity and Eigenvalues

### Characteristic Polynomial

The eigenvalues of the flux Jacobian come from the characteristic polynomial:
```math
P_{2n+1}(x) = Q_n(x) R_{n+1}(x)
```

where $R_{n+1}$ is defined by the HyQMOM closure parameters.

### Eigenvalue Structure

For the 35-moment system, eigenvalues come from:

**1D HyQMOM polynomials**: $Q_2(x)$ and $R_3(x)$ in each direction

**Cross-moment blocks**: Additional eigenvalues from 2D/3D coupling

### Wave Speeds

**Low Mach number**: Eigenvalues ≈ roots of $R_3(x)$

**High Mach number**: Additional fast waves from cross-moment coupling

**Maximum wave speed**: Determines CFL stability condition

## Numerical Implementation

### Moment System Structure

The 35 moments are organized as:
```
[m₀₀₀, m₁₀₀, m₂₀₀, m₃₀₀, m₄₀₀,  # Pure x-direction
 m₀₁₀, m₀₂₀, m₀₃₀, m₀₄₀,        # Pure y-direction  
 m₀₀₁, m₀₀₂, m₀₀₃, m₀₀₄,        # Pure z-direction
 m₁₁₀, m₁₀₁, m₀₁₁,              # 2nd-order cross
 m₂₁₀, m₁₂₀, m₂₀₁, m₁₀₂, m₀₂₁, m₀₁₂,  # 3rd-order cross
 m₂₂₀, m₂₀₂, m₀₂₂, m₃₁₀, m₁₃₀, m₃₀₁, m₁₀₃, m₀₃₁, m₀₁₃,  # 4th-order cross
 m₂₁₁, m₁₂₁, m₁₁₂]              # 4th-order triple cross
```

### Flux Computation

**Convert to standardized moments**

**Apply HyQMOM closures** for 5th-order moments

**Check realizability** and correct if needed

**Compute hyperbolic fluxes** using HLL solver

**Update moments** and repeat realizability check

### Collision Operator

The BGK collision operator relaxes moments toward Maxwellian:
```math
\frac{\partial \mathbf{M}}{\partial t}\bigg|_{\text{collision}} = \frac{\mathbf{G} - \mathbf{M}}{\tau_c}
```

where $\mathbf{G}$ are Maxwellian moments and $\tau_c = \text{Kn}/(2\rho\sqrt{\Theta})$.

## Computational Complexity

### Moment Operations
**Standardization**: $O(N_m)$ where $N_m = 35$ is number of moments

**Realizability check**: $O(N_m^3)$ for matrix determinants

**HyQMOM closure**: $O(N_m)$ polynomial evaluations

### Spatial Discretization
**Grid points**: $N_x \times N_y \times N_z$

**Total unknowns**: $35 \times N_x \times N_y \times N_z$

**Flux computation**: $O(N_m^2)$ per grid point for eigenvalue calculation

### Parallel Scaling
**Domain decomposition**: xy-plane partitioning

**Communication**: Halo exchange for neighboring processors

**Load balancing**: Equal grid points per processor
