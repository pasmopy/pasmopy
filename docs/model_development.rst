Model development
=================

This section walk you through the creation of a mechanistic model from text.

How to use
----------
**Text2Model** is a useful tool to build an ordinary differential equation (ODE) model from a text file describing biochemical reactions.

The text file you need to prepare can be divided into three parts:

#. :ref:`Reaction Layer`
#. :ref:`Observable layer`
#. :ref:`Simulation layer`

.. _Reaction Layer:

Reaction layer
^^^^^^^^^^^^^^

   **description** | **parameters** | **initial conditions**

In the reaction layer, you need to describe biochemical reactions.
Each reaction described in the line number *i* will be converted into *i*\ :sup:`th`\  rate equation.
To specify parameters or initial conditions, you can put those information by using ``|``.

* If you don't specify parameters/initial_conditions, they are initialized to 1 and 0, respectively, and the parameter values will be estimated from experimental data.
* If you want to set a parameter value to 1 and don't want to estimate, you can add ``const`` prefix:

.. code-block:: python

   # The Hill coefficient is fixed to 1.
   TF transcribes mRNA | const n=1

* You can impose parameter constraints by specifying line number in the **parameter** section.

.. code-block:: python
   :linenos:

   # Nucleocytoplasmic Shuttling of DUSP
   DUSP translocates from the cytoplasm to the nucleus
   pDUSP translocates from the cytoplasm to the nucleus |2|

In the example above, you can assume that import and export rates were identical for DUSP (line 2) and pDUSP (line 3).

* If the amount of a model species should be held fixed (never consumed) during simulation, you can add ``fixed`` prefix:

.. code-block:: python

   # [Ligand] will be held fixed to 10.0 during simulation
   Ligand binds Receptor --> LR | kf = 1e-6, kr = 1e-1 | fixed Ligand = 10.0


The available rules can be found at :doc:`modules/reaction_rules`.

You can also supply your own terminology in a reaction rule via:

.. code-block:: python

   from pasmopy import Text2Model

   # Supply "releases" to the reaction rule: "is_dissociated"
   mm_kinetics = Text2Model("michaelis_menten.txt")
   mm_kinetics.register_word("is_dissociated", "releases")
   # Now you can use "releases" in your text, e.g., 'ES releases E and P'
   mm_kinetics.convert()

.. _Observable Layer:

Observable layer (Prefix: ``@obs``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the observable layer, you need to specify ``biomass.observables``, which can correlate model simulations and experimental measurements.
You can create an observable by using model parameters (``p``) and species (``u``).
For example, if the total amount of SOS bound to EGFR should be the sum of RGS (EGFR-Grb2-SOS) and RShGS (EGFR-Shc-Grb2-SOS) complexes in your model, then you can write as follows:

.. code-block:: python

   @obs Total_SOS_bound_to_EGFR: u[RGS] + u[RShGS]

.. _Simulation Layer:

Simulation layer (Prefix: ``@sim``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the simulation layer, you can set simulation conditions, e.g, the simulation time span, the initial concentration of model species, etc.

Example:

.. code-block:: python

   @sim tspan: [0, 120]
   @sim unperturbed: init[EGF] = 0
   @sim condition EGF20nM: init[EGF] = 680
   @sim condition EGF2nM: init[EGF] = 68

* **tspan**:

   Two element vector ``[t0, tf]`` specifying the initial and final times.

* **unperturbed**:

   Description of the untreated condition to find the steady state.

* **condition**:

   Experimental conditions. Use ``p`` and ``init`` to modify model parameters and initial conditions, respectively.


Examples
--------

Michaelis-Menten enzyme kinetics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example shows you how to build a simple Michaelis-Menten two-step enzyme catalysis model with Pasmopy.

   E + S ⇄ ES → E + P

An enzyme, E, binding to a substrate, S, to form a complex, ES, which in turn releases a product, P, regenerating the original enzyme.

#. Prepare a text file describing biochemical reactions (``michaelis_menten.txt``)
   
   .. code-block:: python
      :linenos:

      E binds S --> ES | kf=0.003, kr=0.001 | E=100, S=50
      ES dissociates to E and P | kf=0.002, kr=0

      @obs Substrate: u[S]
      @obs E_free: u[E]
      @obs E_total: u[E] + u[ES]
      @obs Product: u[P]
      @obs Complex: u[ES]

      @sim tspan: [0, 100]

#. Convert text into an executable model

   .. code-block:: python

      from pasmopy import Text2Model

      description = Text2Model("michaelis_menten.txt")
      description.convert()

#. Run simulation with biomass_

   .. code-block:: python

      from biomass import Model, run_simulation
      import michaelis_menten

      model = Model(michaelis_menten.__package__).create()
      run_simulation(model)

.. image:: https://raw.githubusercontent.com/pasmopy/pasmopy/master/docs/_static/img/michaelis_menten_sim.png

EGF signaling
^^^^^^^^^^^^^
Below is an example of Pasmopy in action to illustrate EGF signalling pathway. 

Reference:

   Kholodenko, B. N., Demin, O. V, Moehren, G. & Hoek, J. B. Quantification of short term signaling by the epidermal growth factor receptor. *J. Biol. Chem.* **274**, 30169–30181 (1999). https://doi.org/10.1074/jbc.274.42.30169

#. Prepare a text describing EGF signaling in hepatocytes (``Kholodenko_JBC_1999.txt``)

   .. code-block:: python
      :linenos:
      
      EGF binds EGFR --> Ra | kf=0.003, kr=0.06 | EGFR=100
      Ra dimerizes --> R2 | kf=0.01, kr=0.1
      R2 is phosphorylated --> pR2 | kf=1, kr=0.01
      pR2 is dephosphorylated --> R2 | V=450, K=50
      pR2 binds PLCg --> RPL | kf=0.06, kr=0.2 | PLCg=105
      RPL is phosphorylated --> pRPL | kf=1, kr=0.05
      pRPL is dissociated into pR2 and pPLCg | kf=0.3, kr=0.006
      pPLCg is dephosphorylated --> PLCg | V=1, K=100
      pR2 binds Grb2 --> RG | kf=0.003, kr=0.05 | Grb2=85
      RG binds SOS --> RGS | kf=0.01, kr=0.06 | SOS=34
      RGS is dissociated into pR2 and GS | kf=0.03, kr=4.5e-3
      GS is dissociated into Grb2 and SOS | kf=1.5e-3, kr=1e-4
      pR2 binds Shc --> RSh | kf=0.09, kr=0.6 | Shc=150
      RSh is phosphorylated --> pRSh | kf=6, kr=0.06
      pRSh is dissociated into pShc and pR2 | kf=0.3, kr=9e-4
      pShc is dephosphorylated --> Shc | V=1.7, K=340
      pRSh binds Grb2 --> RShG | kf=0.003, kr=0.1
      RShG is dissociated into pR2 and ShG | kf=0.3, kr=9e-4
      RShG binds SOS --> RShGS | kf=0.01, kr=2.14e-2
      RShGS is dissociated into ShGS and pR2 | kf=0.12, kr=2.4e-4
      pShc binds Grb2 --> ShG | kf=0.003, kr=0.1
      ShG binds SOS --> ShGS | kf=0.03, kr=0.064
      ShGS is dissociated into pShc and GS | kf=0.1, kr=0.021
      pRSh binds GS --> RShGS | kf=0.009, kr=4.29e-2
      pPLCg is translocated to cytoskeletal or membrane structures --> pPLCg_I | kf=1, kr=0.03

      # observable layer
      @obs Total_phosphorylated_Shc: u[pRSh] + u[RShG] + u[RShGS] + u[pShc] + u[ShG] + u[ShGS]
      @obs Total_Grb2_coprecipitated_with_Shc: u[RShG] + u[ShG] + u[RShGS] + u[ShGS]
      @obs Total_phosphorylated_Shc_bound_to_EGFR: u[pRSh] + u[RShG] + u[RShGS]
      @obs Total_Grb2_bound_to_EGFR: u[RG] + u[RGS] + u[RShG] + u[RShGS]
      @obs Total_SOS_bound_to_EGFR: u[RGS] + u[RShGS]
      @obs ShGS_complex: u[ShGS]
      @obs Total_phosphorylated_PLCg: u[pRPL] + u[pPLCg]

      # simulation layer
      @sim tspan: [0, 120]
      @sim condition EGF20nM: init[EGF] = 680
      @sim condition EGF2nM: init[EGF] = 68

#. Convert text into an executable model

   .. code-block:: python

      from pasmopy import Text2Model

      description = Text2Model("Kholodenko_JBC_1999.txt")
      description.convert()
   
#. Run simulation with biomass_
   
   .. code-block:: python

      from biomass import Model, run_simulation
      import Kholodenko_JBC_1999

      model = Model(Kholodenko_JBC_1999.__package__).create()
      run_simulation(model)


.. _biomass: https://github.com/okadalabipr/biomass