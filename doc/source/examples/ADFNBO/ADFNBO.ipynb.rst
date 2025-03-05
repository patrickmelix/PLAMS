Worked Example
--------------

Initialization
~~~~~~~~~~~~~~

.. code:: ipython3

   from scm.plams import Settings, Molecule, init, AMSJob
   from scm.plams.recipes.adfnbo import ADFNBOJob

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/ADFNBO/plams_workdir.002

Define molecule
~~~~~~~~~~~~~~~

.. code:: ipython3

   # mol = Molecule("methane.xyz")
   def get_molecule(input_string):
      job = AMSJob.from_input(input_string)
      return job.molecule[""]


   mol = get_molecule("""
   System
       Atoms
            C      0.000000      0.000000      0.000000
            H      0.631600      0.631600      0.631600
            H      0.631600     -0.631600     -0.631600
            H     -0.631600      0.631600     -0.631600
            H     -0.631600     -0.631600      0.631600
       End
   End
   """)

Create and run job
~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   s = Settings()
   s.input.AMS.Task = "SinglePoint"
   s.input.ADF.basis.type = "DZP"
   s.input.ADF.xc.lda = "SCF VWN"
   s.input.ADF.relativity.level = "scalar"
   s.adfnbo = ["write", "spherical", "fock"]

   j = ADFNBOJob(molecule=mol, settings=s)
   r = j.run()

::

   [05.03|14:48:57] JOB plamsjob STARTED
   [05.03|14:48:57] JOB plamsjob RUNNING
   [05.03|14:48:59] JOB plamsjob FINISHED
   [05.03|14:48:59] JOB plamsjob SUCCESSFUL

Print results
~~~~~~~~~~~~~

.. code:: ipython3

   lines = r.get_output_chunk(begin="NATURAL BOND ORBITALS (Summary):", end="Charge unit", inc_begin=True, inc_end=True)
   for line in lines:
       print(line)

::

    NATURAL BOND ORBITALS (Summary):

                                                        Principal Delocalizations
              NBO                 Occupancy    Energy   (geminal,vicinal,remote)
    ===============================================================================
    Molecular unit  1  (CH4)
    ------ Lewis --------------------------------------
       1. CR ( 1) C  1             2.00000    -9.77904
       2. BD ( 1) C  1- H  2       1.99927    -0.42240
       3. BD ( 1) C  1- H  3       1.99927    -0.42240
       4. BD ( 1) C  1- H  4       1.99927    -0.42240
       5. BD ( 1) C  1- H  5       1.99927    -0.42240
    ------ non-Lewis ----------------------------------
       6. BD*( 1) C  1- H  2       0.00023     0.46676
       7. BD*( 1) C  1- H  3       0.00023     0.46676
       8. BD*( 1) C  1- H  4       0.00023     0.46676
       9. BD*( 1) C  1- H  5       0.00023     0.46676
      10. RY ( 1) C  1             0.00000     4.16091
      11. RY ( 2) C  1             0.00000    20.97822
      12. RY ( 3) C  1             0.00000     0.62722
      13. RY ( 4) C  1             0.00000     0.99946
      14. RY ( 5) C  1             0.00000     0.62722
      15. RY ( 6) C  1             0.00000     2.08484
      16. RY ( 7) C  1             0.00000     2.16053
      17. RY ( 8) C  1             0.00000     2.30880
      18. RY ( 9) C  1             0.00000     1.79930
      19. RY (10) C  1             0.00000     1.79930
      20. RY ( 1) H  2             0.00049     0.37840
      21. RY ( 2) H  2             0.00000     0.95619
      22. RY ( 3) H  2             0.00000     0.95619
      23. RY ( 4) H  2             0.00000     1.52007
      24. RY ( 1) H  3             0.00049     0.37840
      25. RY ( 2) H  3             0.00000     0.95619
      26. RY ( 3) H  3             0.00000     0.95619
      27. RY ( 4) H  3             0.00000     1.52007
      28. RY ( 1) H  4             0.00049     0.37840
      29. RY ( 2) H  4             0.00000     0.95619
      30. RY ( 3) H  4             0.00000     0.95619
      31. RY ( 4) H  4             0.00000     1.52007
      32. RY ( 1) H  5             0.00049     0.37840
      33. RY ( 2) H  5             0.00000     0.95619
      34. RY ( 3) H  5             0.00000     0.95619
      35. RY ( 4) H  5             0.00000     1.52007
             -------------------------------
                    Total Lewis    9.99709  ( 99.9709%)
              Valence non-Lewis    0.00093  (  0.0093%)
              Rydberg non-Lewis    0.00199  (  0.0199%)
             -------------------------------
                  Total unit  1   10.00000  (100.0000%)
                 Charge unit  1    0.00000
