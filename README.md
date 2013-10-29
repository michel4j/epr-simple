EPR-simple
==========

A Simple event-by-event simulation of the EPR experiment

How it works:
------------
The simulation consists of a Source object, generating particle pairs, to be analyzed at 2 Detection stations.

1) The Source, and Particles:
    Simply generates two tuples each with 3 parameters corresponding to the "hidden variables".
    The source a single parameter 'variant' which by default is 1.0. 
    A particle pair is generated as follows:
     'e'  - an angle common to both particles selected randomly each time from the range [0, 2pi)
     'p1' - a random property of the left particle selected randomly from [0, 1)
     'p2' - a random property of the right particle selected randomly from [0, 1)
    
    The left particle is the tuple (e, p1, variant)
    The right particle is the tuple (e, p2, variant)

2) The Detection Stations:
    Two stations exists named 'Alice' and 'Bob'. Alice will measure the left particle, while Bob will measure the right particle.
    The detection proceeds as follows:
    - A random angle 'a' is selected in the range [0, 2pi)
    - A transformed value 'C' is calculated using the three particle properties and the detector setting 'a' as C = cos(variant*(e-a))**2 - p
    - A threshold value is selected randomly in the range [0, 1). If the absolute value of C is greather than the threshold, the particle is detected in which case the output will be sign(C), otherwise the particle is not detected and the output is 0.
    - The setting 'a' and the output is registed

Statistical Analysis:
--------------------    
At the end of the simulation, the results from each station are saved to a separate file and the statistics and coincidences between Alice and Bob are analyzed and plotted.

The angle difference between Alice and Bob's setting is calculated converted and rounded to the nearest degree. For each angle in the range [0, 2pi) converted and rounded to the nearest degree, we collect all instances where that angle difference was observed and calculate the probabilities:
    - P++ : Both Alice and Bob measured +1
    - P-- : Both Alice and Bob measured -1
    - P+- : Alice got +1 and Bob got -1
    - P-+ : Alice got -1 and Bob got +1
    - A+ : Alice got +1
    - B+ : Bob got +1
From These probabilities we can calcualte the Expectation value:
    - E(a,b) = P++ + P-- - P+- - P-+
    
The results are then plotted. A sample plot after 10,000,000 iterations is shown in epr.png
The results for each station are saved in numpy array files Alice.npy and Bob.npy (open with numpy.load)

Notes:
-----
Note the similarity with the Weih's et al data described in http://arxiv.org/pdf/quant-ph/0606122.pdf (Fig 1). At first I was worried that the maximum of the E(a,b) did not reach 1. However, this is exactly what the experiment shows as well (see Fig 1b of the above paper).
