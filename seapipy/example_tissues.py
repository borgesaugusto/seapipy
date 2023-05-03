import numpy as np
import scipy.stats as sps
import seapipy.lattice_class as lattice_class
import seapipy.surface_evolver as surface_evolver


class ExampleTissues:
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
        se_object = surface_evolver.SurfaceEvolver(self.vertices,
                                                        self.edges,
                                                        self.cells,
                                                        self.densities,
                                                        self.volumes,
                                                        polygonal=False)

        return se_object

    def get_initial_densities(self):
        cell_volumes = self.lattice.get_normally_distributed_volumes(self.cells,
                                                                     means=(self.parameters['cell_v_mean'],),
                                                                     stds=(self.parameters['cell_v_std'],),
                                                                     weights=None)

        densities = self.lattice.get_normally_distributed_densities(self.edges,
                                                                    1,
                                                                    0.1)

        return cell_volumes, densities

    def evolve(self):
        self.se_object.initial_relaxing(evolve_step=10000)
        self.se_object.evolve_relaxing(10, 2500)
        self.se_object.add_vertex_averaging(100)
        self.se_object.change_scale(0.005)
        self.se_object.evolve_relaxing(5, 5000)
        self.se_object.save_one_step(self.parameters['save_dir'], self.parameters['file_name'])
        self.se_object.evolve_relaxing(5, 5000)

    def save_many_steps(self, max_steps: int = 50, step: int = 50):
        self.se_object.save_many_steps(self.parameters['save_dir'], self.parameters['file_name'], max_steps, step)


class NormalFurrow(ExampleTissues):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.new_densities = self.get_new_densities(axis=self.parameters["axis"])
        self.se_object.change_line_tensions(self.new_densities)
        self.se_object.save_one_step(self.parameters['save_dir'], self.parameters['file_name'])

    def get_new_densities(self, axis="x"):
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

        return new_densities


class CircularFurrow(ExampleTissues):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.volumes, self.densities = self.get_initial_densities()
        self.se_object = self.create_se_file()

        self.se_object.generate_fe_file()
        self.evolve()
        self.new_densities = self.get_new_densities(axis=self.parameters["axis"])
        self.se_object.change_line_tensions(self.new_densities)
        self.se_object.save_one_step(self.parameters['save_dir'], self.parameters['file_name'])

    def get_new_densities(self, axis="x"):
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
        new_densities = {}
        for k, _ in self.edges.items():
            chosen = np.random.choice(self.parameters["edge_tensions"])
            scale = self.parameters["edge_tensions_std"] * chosen

            new_densities[k] = np.random.normal(loc=chosen, scale=scale)

        return new_densities
