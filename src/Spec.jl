using FFTW

@doc raw"""
    Getfreq(dt,Nt)

Returns the frequency vector corresponding to the fft results

```Infos
    dt is the time sample in s
    Nt is the number of samples in the trace given to fft()
    The frequency increment is 1/(Nt*dt)
```

# Examples
```julia-repl
julia>  Getfreq(1e-3,2000)
2000-element Vector{Float64}:
    0.5
    1.0
    1.5
    2.0
    2.5
    3.0
    3.5
    4.0
    ⋮
  996.5
  997.0
  997.5
  998.0
  998.5
  999.0
  999.5
 1000.0
```
"""
function Getfreq(dt,Nt)
    f=[1/(Nt*dt)*i for i=1:Nt]
    return f
end

@doc raw"""
    Spec(trace,t0,dt,Nt,tp,ts,\Deltaph)

Returns the amplitude spectrums of the P ans S phases where amplitudes are significant

```Input
    trace: the trace recorded at a given virtual station, starting at -t0, with time sampling dt, Nt samples
    tp: picked P-wave phase (s)
    ts: picked S-wave phase (s)
    \Deltaph: width of the phases (s) 
```

```Output
    Vec2dP: modulus of the fft coeffs of P phase where larger than max(Vec2dp)/10 
    Vec2dS: modulus of the fft coeffs of S phase where larger than max(Vec2dp)/10
    FrequP: corresponding frequency vector to plot P-wave amp. spectrum
    FrequS: corresponding frequency vector to plot S-wave amp. spectrum
```

# Examples
```julia-repl
julia>  (AmpP,AmpS,fP,fS)=Spec(trace1,0.1,1e-3,2000,0.7,1.25,0.35))
Trace length is 2.0 s using 2000 samples
Trace 1 at offset 3000 m is found at index 30
Trace 2 at offset 5000 m is found at index 50
Trace 1:
phase P is picked at 0.4 s
phase S is picked at 0.75 s
phase width is set to 0.35
Trace 2:
phase P is picked at 0.7 s
phase S is picked at 1.25 s
phase width is set to 0.35
+ Figure with 8 subplots
```
"""
function Spec(trace,t0,dt,Nt,tp,ts,Δph)
    
    # define all indexes for P and S phases
    iP=trunc(Int,(tp+t0)/dt)
    iS=trunc(Int,(ts+t0)/dt)
    idelta=trunc(Int,Δph/dt)+1

    # mute necessary segments for both phases
    traceP=copy(trace)
    traceS=copy(trace)
    traceP[1:iP-1].=0.
    traceS[1:iS-1].=0.
    traceP[iP+idelta:Nt].=0.
    traceS[iS+idelta:Nt].=0.

    println("phase P is picked at $tp s")
    println("phase S is picked at $ts s")
    println("phase width is set to $Δph")
    # Compute Fourier coeffs
    ampP=zeros(Nt)
    ampcP=fft(traceP)
    @. ampP=abs(ampcP)

    ampS=zeros(Nt)
    ampcS=fft(traceS)
    @. ampS=abs(ampcS)
    
    # zoom in non zero coeffs.: P
    ind=findall(x->x>(maximum(ampP)/10),ampP)

    ind_halfP=zeros(Int(trunc(length(ind)/2)))
    ind_halfP=ind[1:Int(trunc(length(ind)/2))]
    
    Vec2dP=ampP[ind_halfP]

    # zoom in non zero coeffs.: S 
    indS=findall(x->x>(maximum(ampS)/10),ampS)

    ind_halfS=zeros(Int(trunc(length(ind)/2)))
    ind_halfS=ind[1:Int(trunc(length(ind)/2))]
    
    Vec2dS=ampS[ind_halfS]

    fs=Getfreq(dt,Nt)

    frequP=zeros(length(ind_halfP))
    frequP=fs[ind_halfP]

    frequS=zeros(length(ind_halfS))
    frequS=fs[ind_halfS]

    return Vec2dP,Vec2dS,frequP,frequS
end

@doc raw"""
    ComputeQgraph(amp1,amp2,freq,index,d1,d2,V)

Returns x and y to plot Q-graph for a given phase using the results of Spec() for two different traces from two stations

```Equation
    Can't write it for now (no LateX in online Editor?)
```

# Examples
```julia-repl
julia>  TBD
```
"""
function ComputeQgraph(amp1,amp2,freq,index,d1,d2,V)

    x=zeros(length(index))
    y=zeros(length(index))


    cpt=0
    for i in index
    cpt += 1
    x[cpt]=freq[i]
    y[cpt]=V*log(amp2[i]/amp1[i])/(π*(d2-d1))
    #    Qp[cpt]=-frequP[i]*π*(d2-d1)/(Vp*log(ampP2[i]/ampP[i]))
    end

    return x,y
end
