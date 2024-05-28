# SeismicQ.jl
*Tout ce que vous avez toujours voulu savoir sur le Q.*
## Package Features
- Compute source functions
## Function Documentation: Sources
```@docs
Ricker
```
## Function Documentation: Rheology
For deviatoric rheology, updates take the form of:
```math  
\tau = \theta_\mathrm{s} \dot\varepsilon + \chi_\mathrm{s} \tau^\mathrm{old}
```
For volumetric rheology, updates take the form of:
```math  
P = P^\mathrm{old} + \theta_\mathrm{b} \nabla V + \chi_\mathrm{b} \nabla V^\mathrm{old}
```
```@docs
θs
χs 
θb
χb
```
## Function Documentation: Treatment
```@docs
Spec
ComputeQgraph
Getfreq
GenAttenuatedRicker
PlotReceiverGather
```