# Sidney-Williams-Code
This code, and its corresponding theory booklet, are designed to explore magnetic mirrors under various conditions. Specifically, the code is already pre-built to deal with
the case of a vertical mirror in a gravity field, with, or without ExB drift causing rotation. This is not the only case which can be modelled, any non-relativistic situation
can be modelled with this code, including any additional velocity independent forces. As well, it can deal with any magnetic field, currently it is built with traditional
coil fields, or the approximate symmetric field, it works as well for coils that are not equidistant, though it would require trivial modification (adding L1, and L2
variables into the B-field function). It also works for any number of non-interacting particles, with any masses, charges, or non-relativistic energies. The results of
simulations match very well with traditional theoretical treatment (which is walked through in the theory booklet), and the plotting program is capable of making animations
for a more cinematic analysis. I was not the sole author of this code, Tal Rubin (currently a graduate student at PPPL) added more functionality to the plotting function,
such that the mirror ratio could be visually derived. And the base code was taken from this github: https://github.com/PyPhy/Python/tree/master/Plasma it has been drastically
modified, and changed to focus exclusively on  magnetic mirrors.

This code is by no means "done". I cannot claim that it works for every possible situation, so here are some ideas of what can be added. 
1) Relativity
2) Particle Interactions
3) Non-Coil Magnetic Fields
As I am giving suggestions, you may assume that I would like this to be a "living" code. Modify it how you need, append the theory booklet, and this README and then add 
a folder to the repository with your modification, that way troubleshooting will be easy, and progress will be obvious. I hope that this code proves useful to you.
-Sidney Williams (September 2021)
