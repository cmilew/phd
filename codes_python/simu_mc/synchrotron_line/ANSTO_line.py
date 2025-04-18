import os
import opengate as gate
import numpy as np
import time
import math
from decimal import Decimal
import click


def create_box_vol(simulation, name, mother, size, translation, material, color):
    """Function creating a box volume with given parameters"""

    box_vol = simulation.add_volume("Box", name)
    box_vol.mother = mother
    box_vol.size = size
    box_vol.translation = translation
    box_vol.material = material
    box_vol.color = color


def create_kill_act_vol(
    simulation,
    kill_act_name,
    kill_act_mother,
    kill_act_size,
    kill_act_translation,
    kill_act_col,
):
    """Function creating a box volume killing any particles touching/crossing it"""

    create_box_vol(
        simulation,
        kill_act_name,
        kill_act_mother,
        kill_act_size,
        kill_act_translation,
        "Vacuum",
        kill_act_col,
    )
    kill_act = simulation.add_actor("KillActor", f"kill_act_{kill_act_name}")
    kill_act.mother = kill_act_name


def create_kill_collim(
    sim, collim_name, mother_vol, z_translation, slit_width, slit_height, slit_depth, m
):
    """Function creating a collimator made of 4 joint blocs (top, bottom left and
    right) of total dimension (20 cm + slit_width) x (20 cm + slit_height) x 1 um
    including a slit at the center of 4 blocs with given dimensions in meters."""

    # top volume of the collimator
    y_top_vol = 0.05 * m + slit_height / 2
    width_top_vol = 0.2 * m + slit_width
    create_kill_act_vol(
        sim,
        f"top_vol_collim_{collim_name}",
        mother_vol,
        [width_top_vol, 0.1 * m, slit_depth],
        [0, y_top_vol, z_translation],
        [1, 0, 0, 1],
    )

    # left volume of the collimator
    x_left_stop_vol = 0.05 * m + slit_width / 2
    height_side_vol = 0.1 * m + slit_height
    y_side_vol = -((height_side_vol / 2) - slit_height / 2)
    create_kill_act_vol(
        sim,
        f"left_vol_collim_{collim_name}",
        mother_vol,
        [0.1 * m, height_side_vol, slit_depth],
        [x_left_stop_vol, y_side_vol, z_translation],
        [1, 0, 0, 1],
    )

    # right volume of the collimator
    create_kill_act_vol(
        sim,
        f"right_vol_collim_{collim_name}",
        mother_vol,
        [0.1 * m, height_side_vol, slit_depth],
        [-x_left_stop_vol, y_side_vol, z_translation],
        [1, 0, 0, 1],
    )

    # bottom volume of the collimator
    y_bot_vol = -(0.05 * m + slit_height / 2)
    create_kill_act_vol(
        sim,
        f"bot_vol_collim_{collim_name}",
        mother_vol,
        [slit_width, 0.1 * m, slit_depth],
        [0, y_bot_vol, z_translation],
        [1, 0, 0, 1],
    )


def create_array_of_elements(n_elements, center_to_center, *, offset_to_center=0):
    """Function creating an array of elements separated by a center_to_center distance
    and centered on 0 by default on the x-axis"""
    max_dist = center_to_center * (n_elements - 1)
    return np.linspace(-max_dist / 2, max_dist / 2, n_elements) + offset_to_center


def create_kill_msc_leaves(
    simulation,
    mother_vol,
    n_leaves,
    m,
    z_translation,
    center_to_center,
    *,
    offset_to_center=0,
):
    # gets positions of each strips center
    center_positions = create_array_of_elements(
        n_leaves, center_to_center, offset_to_center=offset_to_center
    )

    for leave_num, pos in enumerate(center_positions):
        create_kill_act_vol(
            simulation,
            f"leave_{leave_num}",
            mother_vol,
            [150e-6 * m, 3e-3 * m, 8e-3 * m],
            [pos, 0, z_translation],
            [0, 0, 1, 1],
        )


def create_phsp(simulation, phsp_name, phsp_mother, phsp_z_translation, m, phsp_attr):
    """function creating a phsp of 20 cm x 20 cm x 1 um made of vacuum at given
    z_translation"""

    # phsp plane
    create_box_vol(
        simulation,
        phsp_name,
        phsp_mother,
        [0.2 * m, 0.2 * m, 1e-6 * m],
        [0, 0, phsp_z_translation],
        "Vacuum",
        [0, 1, 0, 1],
    )

    # phsp actor
    phsp_actor = simulation.add_actor("PhaseSpaceActor", phsp_name)
    phsp_actor.mother = phsp_name
    phsp_actor.attributes = phsp_attr
    phsp_actor.output = os.path.join(
        os.path.dirname(__file__),
        f"output/{phsp_name}.root",
    )
    phsp_actor.debug = False


def run_simulation(n_part):
    # units
    m = gate.g4_units.m
    keV = gate.g4_units.keV
    deg = gate.g4_units.deg

    # Simulation parameters
    N_PARTICLES = n_part
    N_THREADS = 1
    unit_spec_file = keV
    beam_type = "poly_AlAl"
    energy_type = "mono"  # "mono" or "spectrum"

    if energy_type == "spectrum":
        spec_file = "ANSTO_AlAl"
    if energy_type == "mono":
        beam_energy = 48 * keV
    sleep_time = True  # only for parallel simulation on CC
    physics_list = "G4EmLivermorePolarizedPhysics"
    visu = False

    # Beam dimensions defined at BDA
    # beam_width = 35e-3 * m  # beam width for step phantom setup
    beam_width = 35e-3 * m  # setup step phantom
    beam_height = 1.059e-3 * m

    # Dist between source and collimators
    world_length = 35 * m
    margin_world_edge_source = 0.01 * m
    d_source_col = 33.4 * m
    d_source_msc = 33.62 * m
    d_source_end_line_phsp = 34 * m

    # source radius calculated to emcompass half diagonal of BDA + 10 um margin
    margin = 10e-6 * m
    source_rad = (1 / 2) * np.sqrt(beam_width**2 + beam_height**2) + margin
    beam_div = np.arctan(source_rad / d_source_col) * (180 / math.pi) * deg

    # create the simulation
    sim = gate.Simulation()
    db_path = os.path.join(os.path.dirname(__file__), "data/gate_materials.db")
    sim.volume_manager.add_material_database(db_path)

    # main options
    sim.g4_verbose = False
    sim.visu = visu
    sim.visu_type = "vrml"
    sim.check_volumes_overlap = False
    sim.number_of_threads = N_THREADS
    sim.random_seed = "auto"

    # world size
    sim.world.size = [1 * m, 1 * m, world_length]
    sim.world.material = "Vacuum"
    sim.physics_manager.set_production_cut("world", "all", 10 * m)

    # create target volume for source acceptance angle located just after first collim
    z_translation_col = world_length / 2 - margin_world_edge_source - d_source_col
    create_box_vol(
        sim,
        "target_vol",
        "world",
        [beam_width + margin, beam_height + margin, 1e-6 * m],
        [0, 0, z_translation_col - margin],
        "Vacuum",
        [0, 1, 0, 1],
    )

    # point source definition (point source at wiggler)
    source = sim.add_source("GenericSource", "point_source")
    source.particle = "gamma"
    source.mother = "world"
    source.n = N_PARTICLES / sim.number_of_threads
    source.position.type = "point"

    # source parameters
    source.direction.type = "iso"
    source.direction.theta = [0 * deg, beam_div]
    source.direction.phi = [0 * deg, 360 * deg]
    source.position.translation = [0, 0, world_length / 2 - margin_world_edge_source]
    source.direction.acceptance_angle.volumes = ["target_vol"]
    source.direction.acceptance_angle.intersection_flag = True
    source.direction.acceptance_angle.skip_policy = "ZeroEnergy"
    if energy_type == "mono":
        source.energy.type = "mono"
        source.energy.mono = beam_energy
    else:
        source.energy.type = "spectrum_lines"
        spectrum = np.loadtxt(
            os.path.join(os.path.dirname(__file__), f"data/{spec_file}.txt"),
            skiprows=1,
            delimiter="\t",
        )
        source.energy.spectrum_energy = spectrum[:, 0] * unit_spec_file
        source.energy.spectrum_weight = spectrum[:, 1]

    # first collimator shaping beam and killing all particles touching/crossing it
    create_kill_collim(
        sim, "col", "world", z_translation_col, beam_width, beam_height, 1e-6 * m, m
    )

    # creates MSC sides bloc
    # MSC slit width = 125 slits of 50 um + 124 leaves of 350 um (arbiterary choice to
    # have at least 40 mm of beam width covered)
    msc_slit_width = (125 * 50e-6 + 124 * 350e-6) * m
    z_translation_msc = world_length / 2 - margin_world_edge_source - d_source_msc
    create_kill_collim(
        sim, "msc", "world", z_translation_msc, msc_slit_width, 3e-3 * m, 8e-3 * m, m
    )

    # creates MSC 124 leaves
    create_kill_msc_leaves(sim, "world", 124, m, z_translation_msc, 400e-6 * m)
    z_phsp_end_line = (
        world_length / 2 - margin_world_edge_source - d_source_end_line_phsp
    )

    create_phsp(
        sim,
        f"phsp_ANSTO_{beam_type}_{Decimal(N_PARTICLES):.3E}_events",
        "world",
        z_phsp_end_line,
        m,
        ["KineticEnergy", "Weight", "PrePositionLocal", "PreDirectionLocal"],
    )

    # phys
    sim.physics_manager.physics_list_name = physics_list

    # stat actor
    s = sim.add_actor("SimulationStatisticsActor", "stats")
    s.track_types_flag = True
    s.output = os.path.join(
        os.path.join(os.path.dirname(__file__), "output"), "stats.txt"
    )

    # ---------------------------------------------------------------------
    # start simulation
    if sleep_time:
        rd_time = np.random.randint(0, 120, 1)[0]
        time.sleep(rd_time)
    sim.run(start_new_process=False)

    stats = sim.output.get_actor("stats")
    print(stats)


@click.command()
@click.option("--n_part", default=1, help="Number of particle to simulate")
def main(n_part):
    run_simulation(n_part)


if __name__ == "__main__":
    main()
