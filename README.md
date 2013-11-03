EPR-simple
==========

A Simple event-by-event simulation of the EPR experiment

How it works:
------------
The simulation consists of a Source object, generating particle pairs, to be analyzed at 2 Detection stations. The maths of the model can be summarized as:  

        λ = {e, p, s},
        e ∈ [0..2π) 
        p ∈ [0..1)
        s = {1/2, 1}
        e' = e + 2πs
        n  = 2s
        k = √2
        A(a,λ) = sign(-1ⁿ cos n(a − e)) if [½(|-1ⁿ cos n(a − e)| + 1)]ᵏ > p, 0 otherwise
        B(b,λ) = sign(-1ⁿ cos n(a − e')) if [½(|-1ⁿ cos n(a − e')| + 1)]ᵏ > p, 0 otherwise


1) The Source, and Particles:

Simply generates two tuples each with 3 parameters corresponding to the "hidden variables".
The source has a single parameter `spin (s)` which determine the type of particles produced. For spin 1/2 particles such as electrons, s=1/2 for photons s=1.

A particle pair is generated as follows:  

    `e` - an angle common to both particles selected randomly each time from the range [0, 2pi)
    `p` - a property common to both particles selected randomly from [0, 1)
    
The left particle is the tuple `(e, p, s)`
The right particle is the tuple `(e+2πs, p, s)`

2) The Detection Stations:  

Two stations exist named `Alice` and `Bob`. Alice will measure the left particle, while Bob will measure the right particle.

The detection proceeds as follows:  

    - A random angle `a` is selected in the range [0, 2pi). This is the detector setting.
    - A transformed value `C` is calculated using the particle properties and 
      the detector setting `a` as `C = -1ⁿ cos n(a − e)`. The sign of this value, will 
      ultimately determine which channel the particle will be detected at `+1` or `-1`
    - A threshold value value is then calculate from the `C` as
      `C' = [½(|-1ⁿ cos n(a − e)| + 1)]ᵏ, where k=√2`. This value together with the hidden 
      particle property `p`, will determine if the particle goes through the filter. 
      If `C' > p` the particle goes through. Every particle which goes through the 
      filter is detected by one of the two channels.
    - The setting `a` and the output (`+1`, `-1`, or `0`) are registered locally at each station
      and saved in separate files at the end of the simulation.
      

Statistical Analysis:
--------------------    
The statistics are calculated as follows. The angle difference between Alice and Bob's setting is calculated using the first column of their respective output arrays, converted to degrees and rounded to the nearest degree. For each angle in the range [0, 2pi), we collect all instances where that angle difference was observed. Then we count the number of matches/mismatches between the second columns of Alice and Bobs arrays and calculate the probabilities:  

    - P++ : Both Alice and Bob measured +1
    - P-- : Both Alice and Bob measured -1
    - P+- : Alice got +1 and Bob got -1
    - P-+ : Alice got -1 and Bob got +1
    - A+ : Alice got +1
    - B+ : Bob got +1

From These probabilities we can calcualte the Expectation value:  
    
    - E(a,b) = P++ + P-- - P+- - P-+  
    
A sample plot after 10,000,000 iterations is shown in epr.png

Notes:
-----
Each particle is treated separately from the source to detection in a completely local and realistic manner. In fact the whole simulation can be performed on separate computers.

The model reproduces *almost* exactly the QM correlation for both electrons and photons
and matches experimental data very well. See http://arxiv.org/pdf/quant-ph/0606122.pdf (Fig 1).

Obviously the model violates the CHSH inequality.
