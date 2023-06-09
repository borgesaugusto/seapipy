import numpy as np
import scipy.stats as sps
import seapipy.lattice_class as lattice_class
import seapipy.surface_evolver as surface_evolver


class ExampleTissues:
    """
    Parent class for creating tissues

    :param parameters: Dictionary with all the relevant parameters of the system
    :type parameters: dict
    """
    def __init__(self, parameters):
        self.cells = None
        self.edges = None
        self.vertices = None

        self.parameters = parameters

        self.lattice = self.create_lattice()

        self.volumes, self.densities = self.get_initial_densities()
        self.se_object = self.create_se_file()

        self.se_object.generate_fe_file()

        self.evolve()

    def create_lattice(self):
        """
        Create a lattice with the initial Voronoi tessellation

        :return: Lattice class object with the initial Voronoi tessellation
        :rtype: lattice_class.Lattice
        """
        # TODO Use the lattice_class.create_example_lattice() instead of this function
        lattice = lattice_class.Lattice(self.parameters["n_cells_x"], self.parameters["n_cells_y"])

        lattice.generate_voronoi_tessellation(
            lattice.generate_square_seeds(
                standard_deviation=self.parameters["voronoi_seeds_std"],
                spatial_step=self.parameters["voronoi_seeds_step"]
            )
        )
        self.vertices, self.edges, self.cells = lattice.create_lattice_elements()
        return lattice

    def create_se_file(self):
        """
        Create the initial Surface Evolver clean slate to start writing to it

        :return: Surface Evolver object that will be written to disk
        :rtype: surface_evolver.SurfaceEvolver
        """
        se_object = surface_evolver.SurfaceEvolver(self.vertices,
                                                        self.edges,
                                                        self.cells,
                                                        self.densities,
                                                        self.volumes,
                                                        polygonal=False)

        return se_object

    def get_initial_densities(self):
        """
        Create dictionary of cell volumes and membrane densities to be used in the initial Surface Evolver simulation
        using the parameters provided to the class

        :return: cell volumes and membrane densities
        """
        cell_volumes = self.lattice.get_normally_distributed_volumes(self.cells,
                                                                     means=(self.parameters['cell_v_mean'],),
                                                                     stds=(self.parameters['cell_v_std'],),
                                                                     weights=None)

        densities = self.lattice.get_normally_distributed_densities(self.edges,
                                                                    1,
                                                                    0.1)

        return cell_volumes, densities

    def evolve(self):
        """
        Add the initial evolution steps to the Surface Evolver slate

        :return: None
        """
        self.se_object.initial_relaxing(evolve_step=10000)
        self.se_object.evolve_relaxing(10, 2500)
        self.se_object.add_vertex_averaging(100)
        self.se_object.change_scale(0.005)
        self.se_object.evolve_relaxing(5, 5000)
        self.se_object.save_one_step(self.parameters['save_dir'], self.parameters['file_name'])
        self.se_object.evolve_relaxing(5, 5000)

    def save_many_steps(self, max_steps: int = 50, step: int = 50):
        """
        Add to the Surface Evolver slate saving loop *max_steps* times every *step* steps

        :param max_steps: Number of saves to generate
        :type max_steps: int
        :param step: Evolution steps between savings
        :type step: int
        :return: None
        """
        self.se_object.save_many_steps(self.parameters['save_dir'], self.parameters['file_name'], max_steps, step)


class NormalFurrow(ExampleTissues):
    """
    Child class of :class:`ExampleTissues` to create a tissue with a normal distribution of membrane densities in the
    horizontal or vertical axis

    :param parameters: Dictionary with all the relevant parameters of the system
    :type parameters: dict
    """
    def __init__(self, parameters):
        super().__init__(parameters)

        self.new_densities = self.get_new_densities(axis=self.parameters["axis"])
        self.se_object.change_line_tensions(self.new_densities)
        self.se_object.save_one_step(self.parameters['save_dir'], self.parameters['file_name'])

    def get_new_densities(self, axis="x"):
        """
        Create the densities dictionary with the normal distribution in the horizontal or vertical axis

        :param axis: x or y
        :type axis: str
        :return: Dictionary with the new densities to be assigned to the membranes
        :rtype: dict
        """
        tissue_center = np.mean(self.lattice.get_coordinates(self.vertices), axis=1)
        tissue_min = np.min(self.lattice.get_coordinates(self.vertices), axis=1)
        tissue_max = np.max(self.lattice.get_coordinates(self.vertices), axis=1)
        if axis == "x":
            xrange = np.arange(tissue_min[0], tissue_max[0])
            distribution_peak = np.max(sps.norm.pdf(xrange, loc=tissue_center[0], scale=self.parameters['edge_t_std']))
            norm_value = {k: sps.norm.pdf(self.lattice.get_edge_centroid(k, self.vertices, self.edges)[0],
                                          tissue_center[0],
                                          self.parameters['edge_t_std']) / distribution_peak
                          for k, _ in self.edges.items()}

            new_densities = {k: self.densities[k] + self.parameters['edge_t_mean'] * norm_value[k]
                             for k, edge in self.edges.items()}
        elif axis == "y":
            xrange = np.arange(tissue_min[1], tissue_max[1])
            distribution_peak = np.max(sps.norm.pdf(xrange, loc=tissue_center[1], scale=self.parameters['edge_t_std']))
            norm_value = {k: sps.norm.pdf(self.lattice.get_edge_centroid(k, self.vertices, self.edges)[1],
                                          tissue_center[1],
                                          self.parameters['edge_t_std']) / distribution_peak
                          for k, _ in self.edges.items()}
            new_densities = {k: self.densities[k] + self.parameters['edge_t_mean'] * norm_value[k]
                             for k, edge in self.edges.items()}
        else:
            raise NotImplementedError

        return new_densities


class CircularFurrow(ExampleTissues):
    """
    Child class of :class:`ExampleTissues` to create a tissue with a normal distribution of membrane densities at a fixed distance from the center
    of the tissue. Effectively, gives a circular furrow.

    :param parameters: Dictionary with all the relevant parameters of the system
    :type parameters: dict
    """
    def __init__(self, parameters):
        super().__init__(parameters)

        self.volumes, self.densities = self.get_initial_densities()
        self.se_object = self.create_se_file()

        self.se_object.generate_fe_file()
        self.evolve()
        self.new_densities = self.get_new_densities()
        self.se_object.change_line_tensions(self.new_densities)
        self.se_object.save_one_step(self.parameters['save_dir'], self.parameters['file_name'])

    def get_new_densities(self):
        """
        Create the densities dictionary with the normal distribution in a circular furrow

        :return: Dictionary with the new densities to be assigned to the membranes
        :rtype: dict
        """
        new_densities = {}
        tissue_center = np.mean(self.lattice.get_coordinates(self.vertices), axis=1)

        for k, edge in self.edges.items():
            centroid = self.lattice.get_edge_centroid(k, self.vertices, self.edges)
            distance = np.linalg.norm(tissue_center - centroid)
            distribution_peak = sps.norm.pdf(0, loc=0, scale=self.parameters['edge_t_std'])

            new_densities[k] = sps.norm.pdf(distance,
                                            loc=0,
                                            scale=self.parameters['edge_t_std']) / distribution_peak
        # make densities have an average of 1:
        mean_value = np.mean(list(new_densities.values()))
        new_densities = {k: value / mean_value for k, value in new_densities.items()}
        return new_densities


class RandomCellTypes(ExampleTissues):
    """
    Child class to create a tissue with a random assignment of density to each membrane

    :param parameters: Dictionary with all the relevant parameters of the system
    :type parameters: dict
    """
    def __init__(self, parameters):
        super().__init__(parameters)

        self.volumes, self.densities = self.get_initial_densities()
        self.se_object = self.create_se_file()

        self.se_object.generate_fe_file()
        self.evolve()
        self.new_densities = self.get_new_densities()

        self.se_object.change_line_tensions(self.new_densities)
        self.se_object.save_one_step(self.parameters['save_dir'], self.parameters['file_name'])

    def get_new_densities(self):
        """
        Create the densities dictionary using a random assignment of values. It can be controlled using the parameters
        dictionary.

        :return: Dictionary with the new densities to be assigned to the membranes
        :rtype: dict
        """
        new_densities = {}
        for k, _ in self.edges.items():
            chosen = np.random.choice(self.parameters["edge_tensions"])
            scale = self.parameters["edge_tensions_std"] * chosen

            new_densities[k] = np.random.normal(loc=chosen, scale=scale)

        return new_densities
