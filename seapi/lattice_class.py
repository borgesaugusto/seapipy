from dataclasses import dataclass
import numpy as np
import scipy.spatial as spatial
from copy import deepcopy


@dataclass
class Lattice:
    number_cells_x: int
    number_cells_y: int

    def __post_init__(self):
        self.tesselation = None

    def generate_square_seeds(self, standard_deviation=0, spatial_step=1):
        grid_values = ([[(i + np.random.normal(0, standard_deviation)) * spatial_step,
                         (j + np.random.normal(0, standard_deviation)) * spatial_step]
                        for j in range(0, self.number_cells_y) for i in range(0, self.number_cells_x)])
        return grid_values

    def generate_voronoi_tesselation(self, seed_values):
        self.tesselation = spatial.Voronoi(list(seed_values))

        return self.tesselation

    def create_lattice_elements(self):
        new_vertices = {}
        new_cells = {}
        new_edges = {}
        regions = deepcopy(self.tesselation.regions)
        big_edge = []

        cnum = 1

        regions = self.remove_infinite_regions(regions)

        for c in regions:
            temp_for_cell = []
            temp_vertex_for_cell = []
            if len(c) != 0 and -1 not in c:
                # add first to close the cell
                c.append(c[0])
                for ii in range(0, len(c) - 1):
                    temp_big_edge = []
                    x_coordinate = np.around(np.linspace(round(self.tesselation.vertices[c[ii]][0], 3),
                                                         round(self.tesselation.vertices[c[ii + 1]][0], 3), 2), 3)
                    y_coordinate = np.around(self.line_eq(self.tesselation.vertices[c[ii]],
                                                          self.tesselation.vertices[c[ii + 1]],
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
                # calculate cell centroid
                cnum += 1

        return new_vertices, new_edges, new_cells

    def get_middle_cells(self, cells):
        # TODO implement function
        return None

        middle_cells = []
        other_cells = []
        for c in cells:
            cnum = cells.index(c)
            xcm = np.mean([self.tesselation.vertices[ii][0] for ii in c])
            ycm = np.mean([self.tesselation.vertices[ii][1] for ii in c])
            # print(xcm)
            # if xcm > 80:
            #     right_cells.append(-1 * cnum * area_sign)
            # else:
            #     left_cells.append(-1 * cnum * area_sign)
            # take the middle cells
            if xcm > 60 and xcm < 110:
                middle_cells.append(cnum)
            else:
                other_cells.append(cnum)

        return middle_cells, other_cells

    def remove_infinite_regions(self, regions, max_distance=50):
        # if a region has a vertex at more than 50, remove it
        to_delete = []
        for c in regions:
            distances = []
            if len(c) != 0 and -1 not in c:
                for ii in range(0, len(c) - 1):
                    distances.append(np.linalg.norm(self.tesselation.vertices[c[ii]] - self.tesselation.vertices[c[ii + 1]]))
                distances = np.array(distances)
                if np.any(np.where(distances > max_distance, True, False)):
                    to_delete.append(c)

        for c in to_delete:
            regions.remove(c)
        return regions

    def get_cell_area_sign(self, cell, new_vertices):
        return int(np.sign(self.get_cell_area(cell, new_vertices)))

    @staticmethod
    def get_coordinates(vertices):
        x_coordinates = [x[0] for x in vertices.values()]
        y_coordinates = [y[0] for y in vertices.values()]

        return x_coordinates, y_coordinates

    @staticmethod
    def get_edge_centroid(edge_id, new_vertices, new_edges):
        x_coordinates = [new_vertices[ii][0] for ii in new_edges[edge_id]]
        y_coordinates = [new_vertices[ii][1] for ii in new_edges[edge_id]]
        xm = np.mean(x_coordinates)
        ym = np.mean(y_coordinates)
        return xm, ym

    @staticmethod
    def get_normally_distributed_densities(edges, center=1, standard_deviation=0.01):
        return {k: np.random.normal(center, standard_deviation) for k, _ in edges.items()}

    @staticmethod
    def get_vertex_number(vertex, new_vertices):
        if vertex in new_vertices.values():
            vertex_number = list(new_vertices.keys())[list(new_vertices.values()).index(vertex)]
        else:
            if len(new_vertices) > 0:
                vertex_number = max(new_vertices.keys()) + 1
            else:
                vertex_number = 1
            new_vertices[vertex_number] = vertex
        return vertex_number

    @staticmethod
    def get_enum(edge, new_edges):
        if edge in new_edges.values():
            enum = list(new_edges.keys())[list(new_edges.values()).index(edge)]
        elif edge[::-1] in new_edges.values():
            enum = - list(new_edges.keys())[list(new_edges.values()).index(edge[::-1])]
        else:
            if len(new_edges) > 0:
                enum = max(new_edges.keys()) + 1
            else:
                enum = 1
            new_edges[enum] = [edge[0], edge[1]]
        return enum

    @staticmethod
    def line_eq(p0, p1, x):
        """
        returns the value for the y(x) where y is the equation of a line that goes through p0 and p1
        """
        p0 = np.around(p0, 3)
        p1 = np.around(p1, 3)
        m = (p1[1] - p0[1]) / (p1[0] - p0[0])
        return p0[1] + m * (x - p0[0])

    @staticmethod
    def get_cell_area(cell, new_vertices):
        x = [new_vertices[i][0] for i in cell]
        y = [new_vertices[i][1] for i in cell]
        return 0.5 * (np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

