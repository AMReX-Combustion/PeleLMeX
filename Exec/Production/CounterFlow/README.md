## Counterflow Diffusion Flame

Sample test case for running 2-D counterflow diffusion flame. The initial solution is obtained by first
running a "coolflow" case using the input file ``input_coolflow.2d-regt`` without grid refinement to reduce
computational cost (Note ``prob.do_ignition`` is set to 0). The domain will fill with oxidizer and fuel until
a steady-state is reached when the oxidizer and fuel meet at approx. the center of the domain (Visualization
software such as Paraview can be used to determine when steady-state is reached).

The counterflow flame is then ignited using an ignition kernel by setting ``prob.do_ignition = 1`` and 
``peleLM.initDataPlt_patch_flow_variables = true`` in the input file ``input.2d-regt``, which restores 
the steady-state coolflow plotfile specified by ``amr.initDataPlt``. The ignition kernel parameters 
``prob.ignition_SphT`` and ``prob.ignition_SphRad`` are used in ``PeleLMeX_PatchFlowVariables.cpp`` to 
specify the patched kernel temperature and radius, respectively. 