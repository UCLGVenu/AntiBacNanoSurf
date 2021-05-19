# AntiBacNanoSurf
Antibacterial Nanostructured Surfaces Code

Note - all code is Python. Requires numpy and matplotlib.pyplot to run.

Combined Folder contains the codes for the Figure of Merit calculations, combining absorptions and bending energy. There is a readme file in that folder that briefly details the role of each file and the dependencies. 

MaxwellGarnett contains files that carry out the MG approximation on structures. Parameters can be set within the code. GlassVsSi is set up to test different substrates, NoCoat allows 'normalised' absorptions to be found and TesterBase is the simple base loop that can calculate absorptions over a range of paramters. 

Deflection contains a set of slightly depreciated code - better versions are found as part of the Combined folder. These codes extract Young's modulus for an Si pillar of given thickness, and then use a rule-of-mixtures approach to calculate coated moduli. From this, the deflection and energy under a given load are calculated.
