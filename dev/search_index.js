var documenterSearchIndex = {"docs":
[{"location":"#SeismicQ.jl","page":"Home","title":"SeismicQ.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Tout ce que vous avez toujours voulu savoir sur le Q.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Compute source functions","category":"page"},{"location":"#Function-Documentation:-Sources","page":"Home","title":"Function Documentation: Sources","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Ricker","category":"page"},{"location":"#SeismicQ.Ricker","page":"Home","title":"SeismicQ.Ricker","text":"Ricker(t, t₀, 𝑓₀)\n\nCompute the Ricker function for t, t₀ and 𝑓₀.\n\nf = (10 - 20 (pi 𝑓₀ (t - t₀))^2) exp(-π^2 𝑓₀^2 (t - t₀)^2)\n\nExamples\n\njulia>  Ricker(0.0, 0.0, 0.0)\n1.0\n\n\n\n\n\nRicker(x, x₀, t, t₀, 𝑓₀)\n\nCompute the Ricker function with a dirac.\n\nf = delta(x₀) exp(-((x-x₀)^2)20σ₀^2) (10 - 20 (pi 𝑓₀ (t - t₀))^2) exp(-π^20𝑓₀^20(t - t₀)^2)\n\nExamples\n\njulia> Ricker(2., 2., 0.0, 0.0, 0.0)\n0.1353352832366127\n\n\n\n\n\nRicker(x, x₀, t, t₀, 𝑓₀, σ₀)\n\nCompute the Ricker function with 1D spatial support.\n\nf = exp(-((x-x₀)^2)20σ₀^2) (10 - 20 (pi 𝑓₀ (t - t₀))^2) exp(-π^20𝑓₀^20(t - t₀)^2)\n\nExamples\n\njulia> Ricker(2., 0., 0.0, 0.0, 0.0, 1)\n0.1353352832366127\n\n\n\n\n\nRicker(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀)\n\nCompute the Ricker function with 2D spatial support.\n\nf = exp(-((x-x₀)^2)20σ₀^2) (10 - 20 (pi 𝑓₀ (t - t₀))^2) exp(-π^20𝑓₀^20(t - t₀)^2)\n\nExamples\n\njulia> Ricker(2., 0., 2., 0., 0.0, 0.0, 0.0, 1)\n0.01831563888873418\n\n\n\n\n\n","category":"function"},{"location":"#Function-Documentation:-Rheology","page":"Home","title":"Function Documentation: Rheology","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"f_bulk","category":"page"},{"location":"#SeismicQ.f_bulk","page":"Home","title":"SeismicQ.f_bulk","text":"f_bulk(K)\n\nVolumetric part of the elastic constitutive update:  provides the bulk modulus K for use in \n\nΔp = -K v Δt\n\n\n\n\n\n","category":"function"},{"location":"#Function-Documentation:-Treatment","page":"Home","title":"Function Documentation: Treatment","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Spec\nComputeQgraph\nGetfreq\nGenAttenuatedRicker","category":"page"},{"location":"#SeismicQ.Spec","page":"Home","title":"SeismicQ.Spec","text":"Spec(trace,t0,dt,Nt,tp,ts,\\Deltaph)\n\nReturns the amplitude spectrums of the P ans S phases where amplitudes are significant\n\n    trace: the trace recorded at a given virtual station, starting at -t0, with time sampling dt, Nt samples\n    tp: picked P-wave phase (s)\n    ts: picked S-wave phase (s)\n    Δph: width of the phases (s) \n\n    Vec2dP: modulus of the fft coeffs of P phase where they are larger than max(Vec2dp)/10 \n    Vec2dS: modulus of the fft coeffs of S phase where they are larger than max(Vec2dp)/10\n    FrequP: corresponding frequency vector to plot P-wave amp. spectrum\n    FrequS: corresponding frequency vector to plot S-wave amp. spectrum\n\nExamples\n\njulia>  (AmpP,AmpS,fP,fS)=Spec(trace1,0.1,1e-3,2000,0.7,1.25,0.35))\nTrace length is 2.0 s using 2000 samples\nTrace 1 at offset 3000 m is found at index 30\nTrace 2 at offset 5000 m is found at index 50\nTrace 1:\nphase P is picked at 0.4 s\nphase S is picked at 0.75 s\nphase width is set to 0.35\nTrace 2:\nphase P is picked at 0.7 s\nphase S is picked at 1.25 s\nphase width is set to 0.35\n+ Figure with 8 subplots\n\n\n\n\n\n","category":"function"},{"location":"#SeismicQ.ComputeQgraph","page":"Home","title":"SeismicQ.ComputeQgraph","text":"(x,y)=ComputeQgraph(amp1,amp2,freq,index,d1,d2,V)\n\nReturns x and y to plot Q-graph for a given phase using the results of Spec() for two different traces from two stations\n\n    y=V/(π(d2-d1))*ln(amp2(f)/amp1(f)\n    x=f\n\n    Considering two attenuated waves at two different distances from the source (d1<d2), we obtain the following relation:\n    V/(π(d2-d1))*ln(amp2(f)/amp1(f)=-f/Q+V*n/(\\pi(d2-d1))*ln(d1/d2)\n\n    After plotting (plot(x,y)), a flat curve indicates infinite Q (no attenuation). This is what happens with the \"fake\" data used\n    for the developpment of these routines. In an homogeneous Q medium the graph should be linear. In more complex attenuation media,\n    the graph could be non-linear?\n\n    Velocity is considered constant: the current code does not account for possible dispersion...\n\nExamples\n\njulia>  \n\n\n\n\n\n","category":"function"},{"location":"#SeismicQ.Getfreq","page":"Home","title":"SeismicQ.Getfreq","text":"Getfreq(dt,Nt)\n\nReturns the frequency vector corresponding to the fft results\n\n    dt is the time sample in s\n    Nt is the number of samples in the trace given to fft()\n    The frequency increment is 1/(Nt*dt)\n    Vector goes beyond Nyquist Frequency as fft also does. Only the first half of the vectors are useful for us\n\nExamples\n\njulia>  Getfreq(1e-3,2000)\n2000-element Vector{Float64}:\n    0.5\n    1.0\n    1.5\n    2.0\n    2.5\n    3.0\n    3.5\n    4.0\n    ⋮\n  996.5\n  997.0\n  997.5\n  998.0\n  998.5\n  999.0\n  999.5\n 1000.0\n\n\n\n\n\n","category":"function"},{"location":"#SeismicQ.GenAttenuatedRicker","page":"Home","title":"SeismicQ.GenAttenuatedRicker","text":"time_vec,acc_vec = GenAttenuatedRicker(listₓ,Δt,Nt,Vp,Vs,αp,αs,𝑓₀)  ;\n\nFunction that generates a seismic trace (Ricker wavelet) of P and S waves. This wave is attenuated along time and the function creates a matrix of wave amplitude as a function of time and distance to the source.  The function takes as inputs:     listₓ : vector of geophone distances to the source [m]     Δt : time step of the wave signal [s]      Nt : number of time steps of the wave signal      Vp : P-wave velocity [m/s]      Vs : S-wave velocity [m/s]      αp : attenuation factor for P-wave      αs : attenuation factor for S-wave      𝑓₀ : central frequency of the source [Hz]  \n\nand return:     timevec : vector containing the time steps of the received wave signal     accvec  : matrix of wave acceleration at each geophone position\n\nExamples\n\njulia>  time_vec,acc_vec = GenAttenuatedRicker(0:1000:5000,1e-3,2000,7000,4000,2e-4,4e-4,10.0) \n\n\n\n\n\n\n","category":"function"}]
}
