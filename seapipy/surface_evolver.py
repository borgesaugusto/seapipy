from dataclasses import dataclass
import numpy as np
import io


@dataclass
class SurfaceEvolver:
    vertices: dict
    edges: dict
    cells: dict

    density_values: dict
    volume_values: dict

    def __post_init__(self):
        self.fe_file = io.StringIO()
        # self.density_values = list(map(lambda x: round(x, 3), self.density_values.values))
        self.density_values = {key: round(value, 3) for key, value in self.density_values.items()}

    def generate_fe_file(self):
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

        return self.fe_file

    def add_vertex_averaging(self, how_many=1):
        self.fe_file.write(f"V {how_many}; \n")
        return self.fe_file

    def add_refining_triangulation(self, how_many=1):
        self.fe_file.write(f"r {how_many}; \n")
        return self.fe_file

    def change_scale(self, new_scale):
        self.fe_file.write(f"scale := {new_scale}; \n")
        return self.fe_file

    def evolve_system(self, steps=1):
        self.fe_file.write(f"g {steps}; \n")
        return self.fe_file

    def add_t1_swaps(self, max_size=0.1):
        self.fe_file.write(f"t1_edgeswap edge where length < {max_size}; \n")
        return self.fe_file

    def save_one_step(self, output_directory, file_name):
        self.fe_file.write(f"ff := sprintf '{output_directory}/{file_name}%d.dmp',ii; dump ff; ii+=1; \n")
        return self.fe_file

    def save_many_steps(self, output_directory, file_name, max_steps, time_step=1, averaging=1, max_size=0.1):
        self.fe_file.write('while ii < ' + str(max_steps) + ' do { '
                                                            'g ' + str(time_step) + '; V ' + str(averaging) +
                           '; t1_edgeswap edge where length < 0.1; '
                           'ff := sprintf "' + output_directory + file_name + '%d.dmp",ii;'
                                                                              ' dump ff; ii:=ii+1} \n')
        return self.fe_file

    def save_fe_file(self, file_name):
        with open(f'{file_name}', mode='w') as f:
            print(self.fe_file.getvalue(), file=f)
        return True

    def change_line_tensions(self, new_tensions):
        for eid, tension in new_tensions.items():
            self.fe_file.write(f"set edges density {tension} where original == {abs(eid)}; \n")
        return self.fe_file
