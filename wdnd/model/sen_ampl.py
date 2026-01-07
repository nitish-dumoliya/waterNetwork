from amplpy import AMPL
import copy

ampl = AMPL()

ampl.read("wdn_sipopt.mod")
ampl.readData("wdn.dat")

ampl.setOption("solver", "/usr/local/bin/ipopt_sens")
ampl.setOption("ipopt_options", "run_sens yes sens_boundcheck no")
ampl.setOption("presolve", 0)

# ---- Nominal solve
ampl.solve()


l_var = ampl.getVariable("l").getValues().toDict()

arc_max_dia = {}

if self.data_number == 6:
    fixarcs = set(ampl.getSet("fixarcs").getValues().toList())
else:
    fixarcs = set()

for (i, j, d), val in l_var.items():
    if val > 1e-3:
        if self.data_number == 6:
            if (i, j) in fixarcs or (j, i) in fixarcs:
                continue

        if (i, j) not in arc_max_dia:
            arc_max_dia[(i, j)] = d
        else:
            arc_max_dia[(i, j)] = max(arc_max_dia[(i, j)], d)

