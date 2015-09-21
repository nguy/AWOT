from awot.dropsondes.dropsondes import DropSondes

ds = DropSondes()


def plot_skewt(data, **kwargs):
    ds.plot_skewtlogp(data, **kwargs)


def plot_hodograph(data):
    ds.plot_hodograph(data)


def plot_aux(x_value, y_value, **kwargs):
    ds.plot_aux_graph(x_value, y_value, **kwargs)


def run_thermo_calcs(data):
    ds.run_thermo_calcs(data)
    # ds.dry_lift(data)


def run_shear_calcs(data):
    ds.run_shear_calcs(data)


def plot_params():
    ds.generate_parameter_list()
    ds.plot_thermo_calcs()
    ds.plot_shear_calcs()


def plot_wind_barbs(data):
    ds.plot_wind_barbs(data)
