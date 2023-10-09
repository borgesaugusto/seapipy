[![PyPI version](https://badge.fury.io/py/seapipy.svg)](https://pypi.org/project/seapipy/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/seapipy)](https://pypi.org/project/seapipy/)
[![codecov](https://codecov.io/gh/borgesaugusto/seapipy/graph/badge.svg?token=SJFFTX412I)](https://codecov.io/gh/borgesaugusto/seapipy)
# SeapiPy 
#### A Surface evolver API for python


###  Documentation: 
https://seapipy.readthedocs.io/

---

###  Installation

---

###  Usage
To create a simple tissue 10x10 tissue, you can create a lattice object and initialize the vertices, edges and cells of 
the system. Then, you might create values for the cell volumes directly. You could also create normally distributed
tensions for the edges.
```python
import seapipy as sep
lattice = sep.lattice_class.Lattice(10, 10)
vertices, edges, cells = lattice.create_example_lattice()
volume_values = {k: 500 for k, v in cells.items()}
initial_edges_tensions = lattice.get_normally_distributed_densities(edges)
```

Then, you could create the Surface Evolver object using this variables and then initialize the Surface Evolver slate 
where all the functions will be written into, before saving to disk
```python
se_object = sep.surface_evolver.SurfaceEvolver(vertices, 
                                               edges, 
                                               cells,
                                               initial_edges_tensions, 
                                               volume_values, 
                                               polygonal=False)
se_file = se_object.generate_fe_file()
```

The polygonal=False allows curved edges to exist in the tissue. Now, various Surface Evolver functions might be added 
to the file buffer in the *se_file* variable. For example we could add an initial relaxing for the tissue with
```python
se_object.initial_relaxing()
```

Afterwards, we could add a saving function to create a checkpoint in the Surface Evolver simulation using
```python
se_object.save_one_step("path/to/saving/checkpoint", "step_")
```

Which would save the state of the Surface Evolver simulation at *"path/to/saving/chekcpoint"* with name *"step_"* followed
by the number of times it has been saved. 
Finally you could save the whole Surface Evolver slate into the disk and run it using
```python
se_object.save_fe_file("SurfaceEvolverFile")
sep.command.run_evolver("path/to/SurfaceEvolverFile", "path/to/SurfaceEvolverExecutable")
```
