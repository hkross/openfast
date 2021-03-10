Theory
======

DMST models divide the swept area of the turbine rotor into multiple streamtubes parallel
to the inflow direction, as shown in **figure**. Each streamtube is further divided along
the midline of the rotor, separating the upstream and downstream blade sweeps. Within each
streamtube, turbine dynamics are represented by an actuator disk. Separating the upstream
and downstream blade sweeps allows DMST models to account for the passage of downstream
blades through the wake of upstream blades, improving their performance compared to single
streamtube models.

The DMST model implemented in AeroDyn assumes **insert geometry and operating condition 
assumptions.** The model is based on linear momentum and blade element theories, so it also 
assumes **list linear momentum and blade element theories**. 

The theory is applied by first dividing the rotor swept area into multiple streamtubes in
the lateral direction. If desired, the swept area can also be divided in the vertical
direction, and the solver can be applied to each vertical section separately. The induced
velocities in each streamtube are represented as functions of the upstream and downstream 
induction factors. The upstream induced velocity is given as

.. math::
   V = uV_\infty,

the equilibrium induced velocity between the upstream and downstream streamtubes is given as

.. math::
   V_e = (2u-1)V_\infty,

and the downstream induced velocity is given as 

 .. math::
   V\prime = u\prime V_e = u\prime (2u-1)V_\infty.

and applying linear momentum and blade element theories to each streamtube. This yields
an iterative procedure to solve for the induction factor of each actuator disk.