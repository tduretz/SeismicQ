using SeismicQ,Plots

function main()

# Geophone position [m]
listₓ = 0:100:5000;

# Time domain
Δt  = 1e-3
Nt  = 2000

# Central frequency of the source [Hz]
𝑓₀  = 10.

# Velocities
Vp = 7000 # m/s
Vs = 4000 
αp = 2e-4
αs = 4e-4 # (π * f)/ (Q * V)

#Selected distances from source
d1=3000
d2=5000

TimeVec,AccVec=GenAttenuatedRicker(listₓ,Δt,Nt,Vp,Vs,αp,αs,𝑓₀)
#SeismicMatrix=GenMatrix(Vp,Vs)

dt=Δt

#min and max frequencies to compute Q
fqmin=5.0
fqmax=20.0

Tmax=dt*Nt
time=[dt*i for i=1:Nt]

println("Trace length is $Tmax s using $Nt samples")

ind1=Int(trunc(d1/100))
ind2=Int(trunc(d2/100))

rec1=AccVec[ind1,1:Nt]
rec2=AccVec[ind2,1:Nt]

println("Trace 1 at offset $d1 m is found at index $ind1")
println("Trace 2 at offset $d2 m is found at index $ind2")

ampmin=-0.15
ampmax=0.3

p1=plot(time,rec1,xlabel="time (s)",ylabel="acc. (m.s⁻²)",ylim=(ampmin,ampmax),c=:"blue")
p2=plot(time,rec2,xlabel="time (s)",ylabel="acc. (m.s⁻²)",ylim=(ampmin,ampmax),c=:"red")

# second part rec1:
println("Trace 1:")

ampP,ampS,frequP,frequS=Spec(rec1,0.1,dt,Nt,0.4,0.75,0.35)

ymax=maximum(ampP)

p3=plot(frequP,ampP,title="P",xlabel="f(Hz)",ylabel="Amp",ylim=(0,ymax),c=:"blue")
p4=plot(frequS,ampS,title="S",xlabel="f(Hz)",ylabel="Amp",ylim=(0,ymax),c=:"blue")

#plot(p1,p3,p4,layout=(2,2))

# second part rec2:
println("Trace 2:")

ampP2,ampS2,frequP2,frequS2=Spec(rec2,0.1,dt,Nt,0.7,1.25,0.35)

p5=plot(frequP2,ampP2,title="P2",xlabel="f(Hz)",ylabel="Amp",ylim=(0,ymax),c=:"red")
p6=plot(frequS2,ampS2,title="S2",xlabel="f(Hz)",ylabel="Amp",ylim=(0,ymax),c=:"red")


@show frequP
indf=findall(x->x>fqmin && x<fqmax,frequP)
@show indf

#Qp=zeros(length(indf))

xp,yp=ComputeQgraph(ampP,ampP2,frequP,indf,d1,d2,Vp)

p7=plot(xp,yp,xlabel="f(Hz)",ylim=(-1.0,0.),title="Qgraph for P phase",c=:"green")

xs,ys=ComputeQgraph(ampS,ampS2,frequS,indf,d1,d2,Vs)

@show xs,ys

p8=plot(xs,ys,xlabel="f(Hz)",ylim=(-1.0,0.),title="Qgraph for S phase",c=:"green")


plot(p1,p2,p3,p5,p4,p6,p7,p8,layout=(4,2),legend=false)
#plot(p7,legend=false)

end

main()