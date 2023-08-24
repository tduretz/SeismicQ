using FFTW

function Getfreq(dt,Nt)
    f=[1/(Nt*dt)*i for i=1:Nt]
    return f
end

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

    fs=getfreq(dt,Nt)

    frequP=zeros(length(ind_halfP))
    frequP=fs[ind_halfP]

    frequS=zeros(length(ind_halfS))
    frequS=fs[ind_halfS]

    return Vec2dP,Vec2dS,frequP,frequS
end

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
