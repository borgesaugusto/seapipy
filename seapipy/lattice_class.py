from dataclasses import dataclass
import numpy as np
import scipy.spatial as spatial
from copy import deepcopy


@dataclass
class Lattice:
    """
    Lattice class with all lattice properties

    :param number_cells_x: Number of cells in the x-axis to create
    :type  number_cells_x: int
    :param number_cells_y: Number of cells in the x-axis to create
    :type  number_cells_y: int
    """
    number_cells_x: int
    number_cells_y: int

    tessellation: object = None

    def __post_init__(self):
        ...

    def generate_square_seeds(self, standard_deviation: float = 0, spatial_step: int = 1) -> list:
        """
        Generate seeds for tessellation as a rectangular grid

        :param standard_deviation: Move seeds with a normal with this standard deviation
        :type standard_deviation: float
        :param spatial_step:
        :type spatial_step: int
        :return: List of [x,y] coordinates for the seeds
        :rtype: list
        """
        grid_values = ([[(i + np.random.normal(0, standard_deviation)) * spatial_step,
                         (j + np.random.normal(0, standard_deviation)) * spatial_step]
                        for j in range(0, self.number_cells_y) for i in range(0, self.number_cells_x)])
        return grid_values

    def generate_voronoi_tessellation(self, seed_values: list) -> object:
        """
        Generates the voronoi tessellation from the seed values provided

        :param seed_values: [x, y] list of positions to use as seeds for the voronoi tessellation
        :type seed_values: list
        :return: Voronoi tessellation object generated from the seed values
        :rtype: scipy.spatial.Voronoi()
        """
        self.tessellation = spatial.Voronoi(list(seed_values))

        return self.tessellation

    def create_lattice_elements(self) -> tuple:
        """
        Create the vertices, edges and cells from the scipy.spatial.Voronoi() objects

        :return: Vertices, edges and cell's list
        :rtype: tuple
        """
        new_vertices = {}
        new_cells = {}
        new_edges = {}
        regions = deepcopy(self.tessellation.regions)
        big_edge = []

        cnum = 1

        regions = self.remove_infinite_regions(regions)

        for c in regions:
            temp_for_cell = []
            temp_vertex_for_cell = []
            if len(c) != 0 and -1 not in c:
                # add first to close the cell_vertices
                c.append(c[0])
                for ii in range(0, len(c) - 1):
                    temp_big_edge = []
                    x_coordinate = np.around(np.linspace(round(self.tessellation.vertices[c[ii]][0], 3),
                                                         round(self.tessellation.vertices[c[ii + 1]][0], 3), 2), 3)
                    y_coordinate = np.around(self.line_eq(self.tessellation.vertices[c[ii]],
                                                          self.tessellation.vertices[c[ii + 1]],
                                                          x_coordinate), 3)

                    new_edge_vertices = list(zip(x_coordinate, y_coordinate))
                    # add new edges to the global list
                    for v in range(0, len(new_edge_vertices) - 1):
                        v0 = new_edge_vertices[v]
                        v1 = new_edge_vertices[v + 1]

                        vertex_number_1 = self.get_vertex_number(v0, new_vertices)

                        vertex_number_2 = self.get_vertex_number(v1, new_vertices)

                        enum = self.get_enum([vertex_number_1, vertex_number_2], new_edges)

                        temp_big_edge.append(enum)
                        temp_for_cell.append(enum)
                        temp_vertex_for_cell.append(vertex_number_1)
                        temp_vertex_for_cell.append(vertex_number_2)

                    big_edge.append(temp_big_edge)

                area_sign = self.get_cell_area_sign(temp_vertex_for_cell, new_vertices)
                new_cells[-1 * cnum * area_sign] = temp_for_cell
                # new_cells[cnum] = temp_for_cell
                # calculate cell_vertices centroid
                cnum += 1

        return new_vertices, new_edges, new_cells

    def get_middle_cells(self):
        # TODO implement function
        return None

        # middle_cells = []
        # other_cells = []
        # for c in cells:
        #     cnum = cells.index(c)
        #     xcm = np.mean([self.tessellation.vertices[ii][0] for ii in c])
        #     ycm = np.mean([self.tessellation.vertices[ii][1] for ii in c])
        #     # print(xcm)
        #     # if xcm > 80:
        #     #     right_cells.append(-1 * cnum * area_sign)
        #     # else:
        #     #     left_cells.append(-1 * cnum * area_sign)
        #     # take the middle cells
        #     if xcm > 60 and xcm < 110:
        #         middle_cells.append(cnum)
        #     else:
        #         other_cells.append(cnum)
        #
        # return middle_cells, other_cells

    def remove_infinite_regions(self, regions: list, max_distance: float = 50) -> list:
        """
        Remove regions that have a vertex too far away to be part of the tissue. This solves vertices placed in
        'infinity' by the tessellation

        :param regions: List of the voronoi regions
        :type regions: list
        :param max_distance: Maximum distance of a vertex in a polygon
        :type max_distance: float, optional
        :return: List of regions without the offending ones
        :rtype: list
        """
        to_delete = []
        for c in regions:
            distances = []
            if len(c) != 0 and -1 not in c:
                for ii in range(0, len(c) - 1):
                    distances.append(
                        np.linalg.norm(self.tessellation.vertices[c[ii]] - self.tessellation.vertices[c[ii + 1]]))
                distances = np.array(distances)
                if np.any(np.where(distances > max_distance, True, False)):
                    to_delete.append(c)

        for c in to_delete:
            regions.remove(c)
        return regions

    def get_cell_area_sign(self, cell: list, all_vertices: dict) -> int:
        """
        Get the orientation of the polygon for a cell_vertices using the sign of the area through gauss formula

        :param cell: Vertices of the cell_vertices to determine orientation
        :type cell: list
        :param all_vertices: Dictionary with all the vertices in the lattice
        :type all_vertices: dict
        :return: Sign of the cell_vertices area to determine orientation of the polygon
        :rtype: int
        """
        return int(np.sign(self.get_cell_area(cell, all_vertices)))

    def create_example_lattice(self, voronoi_seeds_std: float = 0.15, voronoi_seeds_step: int = 20) -> tuple[dict, dict, dict]:
        """
        Create a lattice with the initial Voronoi tessellation

        :param voronoi_seeds_std: Noise taken from a normal of mean zero and this value as standard deviation
        :type voronoi_seeds_std: float
        :param voronoi_seeds_step: SSpatial step to deposition of seeds in the tessellation
        :type voronoi_seeds_step: float
        :return: The vertices, edges and cells generated
        :rtype: tuple
        """
        self.generate_voronoi_tessellation(
            self.generate_square_seeds(standard_deviation=voronoi_seeds_std,
                                       spatial_step=voronoi_seeds_step))
        vertices, edges, cells = self.create_lattice_elements()

        return vertices, edges, cells

    @staticmethod
    def get_coordinates(vertices: dict) -> tuple:
        """
        Get the coordinates for all vertices in the system

        :param vertices: List of vertices
        :type vertices: dict
        :return: X and Y coordinates
        :rtype: tuple
        """
        x_coordinates, y_coordinates = zip(*vertices.values())

        return list(x_coordinates), list(y_coordinates)

    @staticmethod
    def get_edge_centroid(edge_id: int, vertices: list, edges: list) -> tuple:
        """
        Get the center of the requested edge

        :param edge_id: id of edge
        :type edge_id: int
        :param vertices: List of vertices
        :type vertices: list
        :param edges: list of edges
        :type edges: list
        :return: x and y coordinates of the edge's centroid
        """
        x_coordinates = [vertices[ii][0] for ii in edges[edge_id]]
        y_coordinates = [vertices[ii][1] for ii in edges[edge_id]]
        xm = np.mean(x_coordinates)
        ym = np.mean(y_coordinates)
        return xm, ym

    @staticmethod
    def get_normally_distributed_densities(edges: dict, center: float = 1, standard_deviation: float = 0.01) -> dict:
        """
        Get a dictionary of tensions that is normally distributed around the center with a certain standard deviation

        :param edges: List of edges in the lattice
        :type edges: dict
        :param center: Center of the distribution of tensions
        :type center: float
        :param standard_deviation: Standard deviation for the distribution of tensions
        :type standard_deviation: float
        :return: Dictionary with the keys as edge ids and the tensions as values
        :rtype: dict
        """
        return {k: np.random.normal(center, standard_deviation) for k, _ in edges.items()}

    @staticmethod
    def get_normally_distributed_volumes(cells: dict, means=(500,), stds=(50,), weights=None) -> dict:
        assert len(means) == len(stds)
        if len(means) > 1:
            assert (weights is None) or (
                        len(weights) == len(means)), f'"weights" and "values" have to be of equal length'
            v = np.random.choice(len(means), len(cells), p=weights)
            return {k: int(np.random.normal(means[v[i]], stds[v[i]])) for i, k in enumerate(cells.keys())}
        return {k: int(np.random.normal(means[0], stds[0])) for k in cells.keys()}

    @staticmethod
    def get_vertex_number(vertex: list, vertices: dict) -> int:
        """
        Get id of the vertex from the vertex position or new possible ID if the vertex is already known

        :param vertex: Vertex position
        :type vertex: list
        :param vertices: Dictionary of the vertices currently identified
        :type vertices: dict
        :return: Vertex id if the vertex exists, if not the next possible id to assign
        :rtype: int
        """
        if vertex in vertices.values():
            vertex_number = list(vertices.keys())[list(vertices.values()).index(vertex)]
        else:
            if len(vertices) > 0:
                vertex_number = max(vertices.keys()) + 1
            else:
                vertex_number = 1
            vertices[vertex_number] = vertex
        return vertex_number

    @staticmethod
    def get_enum(edge: list, edges: dict) -> int:
        """
        Get id of the edge from the vertex position or new possible ID if the vertex is already known

        :param edge: id of the vertices that create the edge
        :type edge: list
        :param edges: Dictionary of the edges currently identified
        :type edges: dict
        :return: edge id if the edge exists, if not the next possible id to assign
        :rtype: int
        """
        if edge in edges.values():
            enum = list(edges.keys())[list(edges.values()).index(edge)]
        elif edge[::-1] in edges.values():
            enum = - list(edges.keys())[list(edges.values()).index(edge[::-1])]
        else:
            if len(edges) > 0:
                enum = max(edges.keys()) + 1
            else:
                enum = 1
            edges[enum] = [edge[0], edge[1]]
        return enum

    @staticmethod
    def line_eq(p0: float, p1: float, x: np.ndarray) -> list:
        """
        Get the value of linear function that goes through p0 and p1 at x

        :param p0: First point to join
        :type p0: float
        :param p1: Second point to join
        :type p0: float
        :param x: Independent variable of the equation
        :type p0: list
        :return: Value of y(x) for a line that goes through p0 and p1
        """
        p0 = np.around(p0, 3)
        p1 = np.around(p1, 3)
        m = (p1[1] - p0[1]) / (p1[0] - p0[0])
        return p0[1] + m * (x - p0[0])

    @staticmethod
    def get_cell_area(cell_vertices: list, vertices: dict) -> float:
        """
        Area of a polygon using the shoelace formula

        :param cell_vertices: List of vertex's id of the polygon
        :param vertices: Dictionary of vertices in the system
        :return: Area of the polygon. Positive number indicates clockwise orientation, negative number
        counterclockwise.
        """
        x = [vertices[i][0] for i in cell_vertices]
        y = [vertices[i][1] for i in cell_vertices]
        return 0.5 * (np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
