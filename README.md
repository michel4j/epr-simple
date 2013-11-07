EPR-simple
==========

A Simple event-by-event simulation of the EPR experiment

How it works:
------------
The simulation consists of a Source object, generating particle pairs, to be analyzed at 2 Detection stations. The maths of the model can be summarized as:  

        λ = {e, p, s},  e ∈ [0..2π), p ∈ [0..π/4), s = {1/2, 1}
        e' = e + 2πs
        A(a,λ) = sign(-1ⁿ cos n(a − e)) if |cos n(a − e)| > √2pᵏ, 0 otherwise
        B(b,λ) = sign(-1ⁿ cos n(b − e')) if |cos n(b − e')| > √2pᵏ, 0 otherwise
        where n = 2s, k=π/2

1) The Source, and Particles:

Simply generates two tuples each with 3 parameters corresponding to the "hidden variables".
The source has a single parameter `spin (s)` which determines the type of particles produced. For spin 1/2 particles such as electrons, s=1/2 for photons s=1.

A particle pair is generated as follows:  

    `e` - an angle common to both particles selected randomly each time from the range [0, 2pi)
    `p` - a property common to both particles selected randomly from [0, π/4)
    
The left particle is the tuple `(e, p, n)`
The right particle is the tuple `(e + 2πs, p, n)`

2) The Detection Stations:  

Two stations exist named `Alice` and `Bob`. Alice will measure the left particle, while Bob will measure the right particle.

The detection at each station proceeds as follows:  

    - A random angle `x` is selected in the range [0, 2π). This is the detector setting.
    - A transformed value `C` is calculated using the particle properties and 
      the detector setting `x` as `C = -1ⁿ cos n(x − e)`. The sign of this value, will 
      ultimately determine which channel the particle will be detected at; `+1` or `-1`
    - The absolute value of `C` together with the particle property `p` will determine 
      if the particle goes through the filter.
      If `|C|` > √2p^π/2 the particle goes through. Every particle which goes through the 
      filter is detected by one of the two channels.
    - The setting `x` and the output (`+1`, `-1`, or `0`) are registered locally at each station
      and saved in separate files at the end of the simulation. Each station is not aware of and 
      uses no information from or about the other station. 
      

Statistical Analysis:
--------------------    
The statistics are calculated as follows. The angle difference between Alice and Bob's setting is calculated using the first column of their respective output arrays, converted to degrees and rounded to the nearest degree. For each angle in the range [0, 2π), we collect all instances where that angle difference was observed. Then we count the number of matches/mismatches between the second columns of Alice and Bobs arrays and calculate the probabilities:  

    - N⁺⁺ : Number of pairs where Both Alice and Bob measured +1
    - N⁻⁻ : ...  Both Alice and Bob measured -1
    - N⁺⁻ : ...  Alice got +1 and Bob got -1
    - N⁻⁺ : ...  Alice got -1 and Bob got +1
    - nA⁺ : ...  Alice got +1
    - nB⁺ : ...  Bob got +1

From these counts, we calcualte the individual probabilities:  

    Pⁱʲ = Nⁱʲ/(N⁺⁺ + N⁻⁻ + N⁺⁻ + N⁻⁺), ij ∈ {++, --, +-, -+}
    
The probability for single sided results Aⁱ and Bⁱ are calculated:  

    pAⁱ = nA/(nA⁺ + nA-), pBⁱ = nB/(nB⁺ + nB-), i∈ {+,-}
    
From these probabilities we can calcualte the Expectation value:  

    E(a,b) = P⁺⁺ + P⁻⁻ - P⁺⁻ - P⁻⁺   

The results are then plotted for every angle pair (a,b) in the range [0, 2π). A sample plot after 50,000,000 iterations is shown in `epr.png`. The output for the Bell-test angles (0, 22.5, 45, 67.5) are shown below:  
    
        <a1b1>: E(  0.0, 22.5), AB=-0.93, QM=-0.92
        <a2d2>: E(  0.0, 67.5), AB=-0.39, QM=-0.38
        <c3b3>: E( 45.0, 22.5), AB=-0.92, QM=-0.92
        <c4d4>: E( 45.0, 67.5), AB=-0.93, QM=-0.92
        
        Same Angle <AB> = -1.00, QM = -1.00
        Oppo Angle <AB> = +1.00, QM = +1.00
        CHSH: < 2.0, MODEL: 2.396, QM: 2.389


Notes:
-----
Each particle is treated separately from the source to detection in a completely local and realistic manner. In fact the whole simulation can be performed on separate computers.

Each detection station behaves exactly the same as the other. Swapping the particles and sending them the opposite way does not change the results.

The model reproduces *almost* exactly the QM correlation for both electrons and photons
and matches experimental data very well. See http://arxiv.org/pdf/quant-ph/0606122.pdf (Fig 1).

Obviously the model violates the CHSH inequality.

To run the simulation yourself, you need at least 4GB of memory but you can run fewer iterations by changing the NUM_ITERATIONS constant. 
On my Computer, it runs ~20,000 particle pairs per second and I can get reasonable statistics from about 1,000,000 photon pairs.

On any Linux system, you need to have `matplotlib` installed for the plotting, as well as `numpy`. Then simply run:  
        
        python epr.py

and wait for it to complete. It will show you a progress bar as well as how many particle pairs it is simulating per second. I haven't tested on Windows but you probably want a "batteries-included" python distribution such as "Python xy" http://code.google.com/p/pythonxy/. 
