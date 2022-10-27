ASE AMSCalculator
-------------------

.. currentmodule:: scm.plams.interfaces.adfsuite.ase_calculator

**Example**: See the :ref:`ASECalculatorExample` example.

An interface between PLAMS and ASE is provided through |AMSCalculator|. A
|Settings| object is provided to set up either an |AMSJob| or an |AMSWorker| if
`amsworker = True`. Note that the settings objectcan be overwritten if a
property is requested from the calculator. For example, if `Gradients no` is
present in the `Properties` block, but `forces` are requested in
`AMSCalculator.calculate`, then `Gradients yes` is used. When using `amsworker
= True`, make sure to use |AMSCalculator| as a context manager, or explicitly
call `AMSCalculator.stop_worker()`. When creating a deepcopy of |AMSCalculator|
the |AMSWorker| is not copied and instead every copy has a reference to the
same |AMSWorker|.

AMSCalculator API
~~~~~~~~~~~~~~~~~

.. autoclass:: AMSCalculator
    :exclude-members: __init__, __weakref__, results_from_ams_results
    :no-private-members:

