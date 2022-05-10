import os
from pathlib import Path

import numpy as np
try:
    import dill as pickle
except ImportError:
    import pickle

from scm.plams import AMSJob

def test_hybrid_engine_input():
    """test :meth:`AMSJob.get_input` for writing sub engines in a hybrid engine block of ams."""
    AMSinput = """system
  Atoms
              C      -1.6447506665       1.4391568332       0.0000000000 
              C      -0.5773632247       0.3290989412      -0.0074515501 
              H      -1.2586274675       2.3064141938       0.5149500419 
              H      -2.5299121499       1.0840525749       0.5067446240 
              H      -1.8942698087       1.7054858887      -1.0164689035 
              H       0.3077982588       0.6842031995      -0.5141961741 
              H      -0.9634864236      -0.5381584194      -0.5224015919 
              H      -0.3278440824       0.0627698856       1.0090173534 
  End
  BondOrders
     2 1 1.0
     3 1 1.0
     4 1 1.0
     5 1 1.0
     6 2 1.0
     7 2 1.0
     8 2 1.0
  End
End

task SinglePoint

Engine hybrid
  energy
    term
      engineid ADF-1
      factor 0.5
      region *
    End
    term
      engineid DFTB-1
      factor 0.5
      region *
    End
  End
  engine ADF ADF-1
  EndEngine
  engine DFTB DFTB-1
  EndEngine

EndEngine

"""
    job = AMSJob.from_input(AMSinput)
    assert job.get_input() == AMSinput

def test_list_block_input():
    """test :meth:`AMSJob.get_input` for writing list-blocks"""
    AMSinput = """system
  Atoms
              C      -1.6447506665       1.4391568332       0.0000000000 
  End
End

task SinglePoint

Engine hybrid
  energy
    term
      engineid ADF-1
      factor 0.5
      region *
    End
    term
      engineid ADF-2
      factor 0.5
      region *
    End
    term
      engineid DFTB-1
      factor 0.5
      region *
    End
  End
  engine ADF ADF-1
     Save TAPE10
     Save TAPE11
     Basis
     PerAtomType
     Symbol H
     End
     PerAtomType
     Symbol O
     End
     End
  EndEngine
  engine ADF ADF-2
  EndEngine
  engine DFTB DFTB-1
  EndEngine

EndEngine

"""
    job = AMSJob.from_input(AMSinput)
    assert job.get_input() == AMSinput

