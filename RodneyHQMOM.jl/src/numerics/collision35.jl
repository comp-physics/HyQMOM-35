"""
    collision35(M, dt, Kn)

Apply elastic BGK collision operator to moments.

Relaxes moments toward Maxwellian equilibrium using BGK model:
dM/dt = (MG - M)/tc

# Arguments
- `M`: 35-element moment vector
- `dt`: Time step
- `Kn`: Knudsen number

# Returns
- `Mout`: Updated moments after collision
"""
function collision35(M, dt, Kn)
    # Extract conserved quantities and compute temperature
    rho = M[1]
    umean = M[2] / rho
    vmean = M[6] / rho
    wmean = M[16] / rho
    
    # Compute temperature from trace of covariance matrix
    C200 = M[3]/rho - umean^2
    C020 = M[10]/rho - vmean^2
    C002 = M[20]/rho - wmean^2
    Theta = (C200 + C020 + C002) / 3
    
    # Collision time scale
    tc = Kn / (rho * sqrt(Theta) * 2)
    
    # Maxwellian equilibrium (isotropic covariance)
    MG = InitializeM4_35(rho, umean, vmean, wmean, Theta, 0.0, 0.0, Theta, 0.0, Theta)
    
    # BGK relaxation: dM/dt = (MG - M)/tc
    Mout = MG - exp(-dt/tc) * (MG - M)
    
    return Mout
end
