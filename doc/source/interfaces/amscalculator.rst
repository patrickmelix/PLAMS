ASE AMSCalculator
-------------------

.. currentmodule:: scm.plams.interfaces.adfsuite.ase_calculator

**Example**: See the :ref:`ASECalculatorExample` example.

An interface between PLAMS and ASE is provided through |AMSCalculator|. A
|Settings| object is used to set up either an |AMSJob| or an |AMSWorker| if
`amsworker = True`. The |Settings| object defines which ASE properties are available:

.. list-table::
   :widths: 25 25
   :header-rows: 1

   * - ASE property
     - AMS settings
   * - `energy`
     - no setting required
   * - `forces`
     - `properties.gradients`
   * - `stress`
     - `properties.stresstensor`

When using `amsworker = True`, make sure to use |AMSCalculator| as a context manager, or explicitly
call `AMSCalculator.stop_worker()`. When creating a deepcopy of |AMSCalculator|
the |AMSWorker| is not copied and instead every copy of |AMSCalculator| has a reference to the
same |AMSWorker|.

AMSCalculator API
~~~~~~~~~~~~~~~~~

.. autoclass:: AMSCalculator
    :exclude-members: __init__, __weakref__, results_from_ams_results
    :no-private-members:

