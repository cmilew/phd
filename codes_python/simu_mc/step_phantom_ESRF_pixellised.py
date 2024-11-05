import os
import opengate as gate
import numpy as np
import time
import uproot
import click

# units
m = gate.g4_units.m
cm = gate.g4_units.cm
mm = gate.g4_units.mm
um = gate.g4_units.um
nm = gate.g4_units.nm
keV = gate.g4_units.keV
deg = gate.g4_units.deg


def create_box_vol(
    simulation,
    box_name,
    box_mother,
    box_size,
    box_translation,
    box_material,
    box_color,
    box_cut,
):
    """Function creating a box volume and setting its properties"""

    box = simulation.add_volume("Box", box_name)
    box.mother = box_mother
    box.size = box_size
    box.translation = box_translation
    box.material = box_material
    box.color = box_color
    simulation.physics_manager.set_production_cut(box_name, "all", box_cut)


def run_simulation(n_part):
    # simulation parameters
    N_PARTICLES = n_part
    N_THREADS = 1
    visu = False
    sleep_time = False
    phsp_filename = "phsp_esrf_line_1.000E+5_events.root"
    physics_list = "G4EmLivermorePolarizedPhysics"

    # setup parameters
    nb_offset_strips = 0
    # negative = towards 136, positive = towards 1
    direction_offset = 1
    slab_material = "RW3"
    slab_w = 16 * cm
    slab_h = 2 * cm
    slab_d = 0 * cm
    cut_slab = 1 * mm
    cut_detector = 10 * um

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
    sim.world.size = [1 * m, 1 * m, 7 * m]
    sim.world.material = "G4_AIR"
    sim.physics_manager.set_production_cut("world", "all", 1 * m)

    # phsp plane
    create_box_vol(
        sim,
        "phsp_plane",
        "world",
        [20 * cm, 10 * cm, 1 * um],
        [0, 0, 3.5 * m],
        "G4_AIR",
        [1, 0, 0, 1],
        1 * mm,
    )

    # phsp source
    source = sim.add_source("PhaseSpaceSource", "phsp_source")
    source.mother = "phsp_plane"
    source.particle = "gamma"
    source.n = N_PARTICLES / sim.number_of_threads
    source.batch_size = 100_000
    source.position.translation = [0, 0, 0 * m]
    source.phsp_file = os.path.join(
        os.path.dirname(__file__), f"output/{phsp_filename}"
    )
    source.position_key = "PrePositionLocal"
    source.direction_key = "PreDirectionLocal"
    source.energy_key = "KineticEnergy"
    source.weight_key = "Weight"
    if sim.number_of_threads > 1:
        source.entry_start = [e * source.n for e in range(sim.number_of_threads)]
        phsp = uproot.open(source.phsp_file)
        first_branch = phsp.keys()[0]
        phsp_n = int(phsp[first_branch].num_entries)

    # slab definition
    if slab_d != 0:
        create_box_vol(
            sim,
            "rw3_slab",
            "world",
            [slab_w, slab_h, slab_d],
            [0, 0, 2.5 * m],
            slab_material,
            [0, 0, 1, 1],
            cut_slab,
        )

    # Diamond detector definition
    # lateral shift of det to define which strip is facing the mb
    lat_shift = (232.5 / 2 + 232.5 * nb_offset_strips) * direction_offset

    # 136 strips * 0.1725 mm + 135 interstrips of 0.06 mm + 2 sides of diamond of 0.06mm
    # = 31.68 mm of detector width
    create_box_vol(
        sim,
        "diamond_detector",
        "world",
        [31.68, 3.32 * mm, 150e-6 * m],
        [lat_shift * um, 0, -2.8 * m],
        "diamant_det",
        [0, 1, 1, 1],
        cut_detector,
    )

    # volume containing all the strips for dose actor retrievements
    # 136 strips of 172.5 um and 135 interstrips of 60 um = 31.56 mm
    create_box_vol(
        sim,
        "dose_actor_vol",
        "diamond_detector",
        [31.56 * mm, 3.2 * mm, 150e-6 * m],
        [0, 0, 0],
        "diamant_det",
        [1, 0, 0, 1],
        cut_detector,
    )

    # dose actor
    dose = sim.add_actor("DoseActor", "dose_detector")
    dose.mother = "dose_actor_vol"
    dose.output = os.path.join(
        os.path.join(os.path.dirname(__file__), "output"),
        f"{slab_d}cm_{slab_material}_detector.mhd",
    )
    dose.square = False
    dose.dose = False
    # pixel width = 0.0075 mm => 23 pixels per strips and 8 pixels per interstrip space
    # 23 * 136 strips + 8 * 135 interstrips = 4208 pixels
    dose.size = [4208, 1, 1]
    dose.spacing = [0.0075 * mm, 3.2 * mm, 150 * um]

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
    print(f"source.n = {source.n}")
    sim.run(start_new_process=False)

    stats = sim.output.get_actor("stats")
    print(stats)


@click.command()
@click.option("--n_part", default=1, help="Number of particle to simulate")
def main(n_part):
    run_simulation(n_part)


if __name__ == "__main__":
    main()
