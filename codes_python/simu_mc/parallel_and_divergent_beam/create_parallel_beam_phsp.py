import os
import opengate as gate
import numpy as np
import time
from decimal import Decimal
import click


def create_box_vol(sim, name, mother, size, translation, material, color):
    """Function creating a box volume with given parameters"""

    box_vol = sim.add_volume("Box", name)
    box_vol.mother = mother
    box_vol.size = size
    box_vol.translation = translation
    box_vol.material = material
    box_vol.color = color


def create_phsp(sim, name, mother, x_dim, y_dim, z_pos, m, attr):
    """function creating a phsp of x_dim * y_dim * 1 um made of vacuum at z_pos"""

    # phsp plane
    create_box_vol(
        sim,
        f"plane_{name}",
        mother,
        [x_dim, y_dim, 1e-6 * m],
        [0, 0, z_pos],
        "Vacuum",
        [1, 0, 0, 1],
    )

    # phsp actor
    phsp_actor = sim.add_actor("PhaseSpaceActor", name)
    phsp_actor.mother = f"plane_{name}"
    phsp_actor.attributes = attr
    phsp_actor.output = os.path.join(
        os.path.dirname(__file__),
        f"output/{name}.root",
    )
    phsp_actor.debug = False


def run_simulation(n_part):
    # units
    m = gate.g4_units.m
    # eV = gate.g4_units.eV
    keV = gate.g4_units.keV

    # Simulation parameters
    N_PARTICLES = n_part
    N_THREADS = 1
    # unit_spec_file = eV
    sleep_time = True  # only for parallel simulation on CC
    physics_list = "G4EmLivermorePolarizedPhysics"
    visu = False
    beam_width = 35e-3 * m  # beam width for step phantom setup
    beam_height = 795e-6 * m
    beam_energy = 106 * keV
    world_length = 0.03 * m

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

    # source parameters
    source = sim.add_source("GenericSource", "point_source")
    source.particle = "gamma"
    source.mother = "world"
    source.n = N_PARTICLES / sim.number_of_threads
    source.position.type = "box"
    source.position.size = [beam_width, beam_height, 0.01 * m]
    source.position.translation = [0, 0, 0]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]
    source.energy.type = "mono"
    source.energy.mono = beam_energy

    # phase space plane
    create_phsp(
        sim,
        f"parallel_beam_phsp_{Decimal(N_PARTICLES):.3E}",
        "world",
        beam_width,
        beam_height,
        0.01 * m,
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
