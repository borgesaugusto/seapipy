from dataclasses import dataclass
import io


@dataclass
class SurfaceEvolver:
    """
    Class to interact with the Surface Evolver file buffer.

    :param vertices: Vertices in the tesselation
    :type vertices: dict
    :param edges: Edges in the tesselation
    :type edges: dict
    :param vertices: Cells in the tesselation
    :type cells: dict
    :param density_values: Initial line tensions for the edges
    :type density_values: dict
    :param volume_values: Initial target areas in the tesselation
    :type volume_values: dict
    :param polygonal: Whether to use polygons or allowed curved edges
    :type polygonal: bool, optional
    """
    vertices: dict
    edges: dict
    cells: dict

    density_values: dict
    volume_values: dict

    polygonal: bool = True

    def __post_init__(self):
        self.fe_file = io.StringIO()
        self.density_values = {key: round(value, 3) for key, value in self.density_values.items()}

    def generate_fe_file(self) -> io.StringIO:
        """
        Generate the initial Surface Evolver slate to write into

        :return: Initialized Surface Evolver slate
        :rtype: io.StringIO()
        """
        self.fe_file.write("SPACE_DIMENSION 2 \n")
        self.fe_file.write("SCALE 0.005 FIXED\n")
        self.fe_file.write("STRING \n")
        self.fe_file.write("\n")
        self.fe_file.write("vertices \n")

        for k, v in self.vertices.items():
            self.fe_file.write(f"{k}   {round(v[0], 3)} {round(v[1], 3)}\n")
        self.fe_file.write("\n")
        self.fe_file.write("edges \n")
        for k, v in self.edges.items():
            # lambda_val = round(new_edges_tensions[[k in ee for ee in bedge].index(True)], 3)
            lambda_val = self.density_values[k]
            self.fe_file.write(f"{abs(k)}   {v[0]}   {v[1]}   density {lambda_val}\n")

        self.fe_file.write("\n")
        self.fe_file.write("faces \n")
        for k, v in self.cells.items():
            # color = np.random.choice(list(color_list))
            str_value = " ".join(str(vv) for vv in v)
            self.fe_file.write(f"{abs(k)}   {str_value} \n")

        self.fe_file.write("\n")
        self.fe_file.write("bodies \n")
        for k, v in self.cells.items():
            self.fe_file.write(f"{abs(k)}   {k}    VOLUME {self.volume_values[k]} \n")

        self.fe_file.write("\n \n")
        self.fe_file.write("read \n \n")
        self.fe_file.write("show_all_edges off \n")
        # f.write("clipped on \n")
        self.fe_file.write("metric_conversion off \n")
        self.fe_file.write("autorecalc on \n")
        self.fe_file.write("gv_binary off \n")
        self.fe_file.write("gravity off \n")
        self.fe_file.write("ii := 0; \n")

        if not self.polygonal:
            self.add_refining_triangulation(3)
        return self.fe_file

    def add_vertex_averaging(self, how_many: int = 1) -> io.StringIO:
        """
        Add vertex averaging using the V Surface Evolver function, at the end of the Surface Evolver slate

        :param how_many: Number of times the averaging should be done
        :return: Current Surface Evolver slate
        :rtype: io.StringIO()
        """
        self.fe_file.write(f"V {how_many}; \n")
        return self.fe_file

    def add_refining_triangulation(self, how_many: int = 1) -> io.StringIO:
        """
        Add a mesh refinement using the r Surface Evolver function, at the end of the Surface Evolver slate

        :param how_many: Number of times the refinement should be done
        :type how_many: int
        :return: Current Surface Evolver slate
        :rtype: io.StringIO()
        """
        self.fe_file.write(f"r {how_many}; \n")
        return self.fe_file

    def change_scale(self, new_scale: float) -> io.StringIO:
        """
        Change the scale of the Surface Evolver simulation at the end of the Surface Evolver slate

        :param new_scale: Numerical value of the new scale
        :type new_scale: float
        :return: Current Surface Evolver slate
        :rtype: io.StringIO()
        """
        self.fe_file.write(f"scale := {new_scale}; \n")
        return self.fe_file

    def evolve_system(self, steps: int = 1) -> io.StringIO:
        """
        Add evolution function using the Surface Evolver go function

        :param steps: Number of steps to evolve for
        :type steps: int
        :return: Current Surface Evolver slate
        :rtype: io.StringIO()
        """
        self.fe_file.write(f"g {steps}; \n")
        return self.fe_file

    def add_t1_swaps(self, max_size: float = 0.1) -> io.StringIO:
        """
        Add check for T1 swaps using the Surface Evolver t1_edgeswap function

        :param max_size: Maximum size of interfaces before making a T1 swap
        :type max_size: float
        :return: Current Surface Evolver slate
        :rtype: io.StringIO()
        """
        self.fe_file.write(f"t1_edgeswap edge where length < {max_size}; \n")
        return self.fe_file

    def initial_relaxing(self, evolve_step: int = 2500, averaging: int = 100) -> io.StringIO:
        """
        Initial standard relaxing with vertex averaging and scale change followed by evolution

        :param evolve_step: Number of steps to evolve at each scale chnage
        :type evolve_step: int
        :param averaging: Number of vertex averagings to perform
        :type evolve_step: int
        :return: Current Surface Evolver slate
        :rtype: io.StringIO
        """
        self.add_vertex_averaging(averaging)
        self.change_scale(0.25)
        self.evolve_system(evolve_step)
        self.add_vertex_averaging(averaging)
        self.change_scale(0.1)
        self.evolve_system(evolve_step)
        self.add_vertex_averaging(averaging)
        self.change_scale(0.01)
        return self.fe_file

    def evolve_relaxing(self, number_of_times: int = 1, steps: int = 1, max_size: float = 0.1) -> io.StringIO:
        """
        Evolve the system a fixed number of steps and perform T1 swaps after a definite number of times

        :param number_of_times: Number of times the evolution plus T1 swaps happen
        :type number_of_times: int
        :param steps: Number of evolution step for all iterations
        :type steps: int
        :param max_size: Maximum size allowed for membranes before a T1 happens
        :type max_size: float
        :return: Current Surface Evolver slate
        :rtype: io.StringIO
        """
        for _ in range(0, number_of_times):
            self.evolve_system(steps)
            self.add_t1_swaps(max_size)
        return self.fe_file

    def save_one_step(self, output_directory: str, file_name: str) -> io.StringIO:
        """
        Add a savepoint to the Surface Evolver slate using the sprintf function

        :param output_directory: Folder to save the file into
        :type output_directory: str
        :param file_name: Name of the file to be saved
        :type file_name: str
        :return: Current Surface Evolver slate
        :rtype: io.StringIO
        """
        self.fe_file.write(f'ff := sprintf "{output_directory}/{file_name}%d.dmp",ii; dump ff; ii+=1; \n')
        return self.fe_file

    def save_many_steps(self, output_directory, file_name, max_steps, time_step=1, averaging=1, max_size=0.1):
        self.fe_file.write('while ii < ' + str(max_steps) + ' do { '
                                                            'g ' + str(time_step) + '; V ' + str(averaging) +
                           '; t1_edgeswap edge where length < 0.1; '
                           'ff := sprintf "' + output_directory + "/" + file_name + '%d.dmp",ii;'
                                                                                    ' dump ff; ii:=ii+1} \n')
        return self.fe_file

    def save_fe_file(self, file_name: str) -> bool:
        """
        Save the Surface Evolver slate to disk

        :param file_name:
        :return: Success state of the saving
        :rtype: bool
        """
        self.fe_file.write('q; \n')
        with open(f'{file_name}', mode='w') as f:
            print(self.fe_file.getvalue(), file=f)
        return True

    def change_line_tensions(self, new_tensions: dict) -> io.StringIO:
        """
        Set new densities for the system's  membranes. Setting is done using the dictionary key as membrane id and the
        value as the new tension

        :param new_tensions:
        :return: Current Surface Evolver slate
        :rtype: io.StringIO
        """
        for eid, tension in new_tensions.items():
            self.fe_file.write(f"set edges density {tension} where original == {abs(eid)}; \n")
        return self.fe_file
