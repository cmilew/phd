import os
import opengate as gate
import numpy as np
import time
import uproot
import click


def create_box_vol(
    simulation,
    box_name,
    box_mother,
    box_size,
    box_translation,
    box_material,
    box_color,
    box_cut_gamma,
    box_cut_electron,
):
    """Function creating a box volume and setting its properties"""

    box = simulation.add_volume("Box", box_name)
    box.mother = box_mother
    box.size = box_size
    box.translation = box_translation
    box.material = box_material
    box.color = box_color
    simulation.physics_manager.set_production_cut(box_name, "gamma", box_cut_gamma)
    simulation.physics_manager.set_production_cut(
        box_name, "electron", box_cut_electron
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
        "G4_AIR",
        [0, 0, 1, 1],
        1e-3 * m,
        1e-3 * m,
    )

    # phsp actor
    phsp_actor = simulation.add_actor("PhaseSpaceActor", phsp_name)
    phsp_actor.attached_to = phsp_name
    phsp_actor.attributes = phsp_attr
    phsp_actor.output_filename = os.path.join(
        os.path.dirname(__file__),
        f"output/{phsp_name}.root",
    )
    phsp_actor.debug = False


def run_simulation(n_part):
    # units
    m = gate.g4_units.m

    # simulation parameters
    N_PARTICLES = n_part
    N_THREADS = 18
    visu = False
    sleep_time = False
    phsp_filename = "phsp_end_ESRF_line1.000E+9_events.root"
    physics_list = "G4EmLivermorePolarizedPhysics"

    world_length = 7 * m
    margin_world_edge_source = 0.01 * m
    # phsp source at 41 m distance from synchrotron source
    d_source_slab = 1 * m
    d_source_det = 6.3 * m

    # Nb of strips offset and direction (negative = towards 136, positive = towards 1)
    nb_offset_strips = 0
    direction_offset = 1

    # slab parameters
    slab_material = "RW3"
    slab_w = 0.16 * m
    slab_h = 0.02 * m
    slab_d = 0 * m

    # cuts
    cut_slab_photon = 1e-3 * m
    cut_detector_photon = 10e-6 * m

    # create the simulation
    sim = gate.Simulation()
    db_path = os.path.join(os.path.dirname(__file__), "../data/gate_materials.db")
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
    sim.world.material = "G4_AIR"
    sim.physics_manager.set_production_cut("world", "all", 1 * m)

    # phsp source plane
    create_box_vol(
        sim,
        "phsp_source_plane",
        "world",
        [0.2 * m, 0.2 * m, 1e-6 * m],
        [0, 0, world_length / 2 - margin_world_edge_source],
        "G4_AIR",
        [1, 0, 0, 1],
        1e-3 * m,
        1e-3 * m,
    )

    # phsp source
    source = sim.add_source("PhaseSpaceSource", "phsp_source")
    source.attached_to = "phsp_source_plane"
    phsp_path = os.path.join(os.path.dirname(__file__), f"../data/{phsp_filename}")
    source.phsp_file = phsp_path
    source.position_key = "PrePositionLocal"
    source.direction_key = "PreDirectionLocal"
    source.energy_key = "KineticEnergy"
    source.weight_key = "Weight"
    source.particle = "gamma"
    source.n = N_PARTICLES / sim.number_of_threads
    source.batch_size = 100_000

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
            [0, 0, world_length / 2 - margin_world_edge_source - d_source_slab],
            slab_material,
            [0, 0, 1, 1],
            cut_slab_photon,
            1 * m,
        )

    # # phsp in front of detector
    # z_phsp_front_msc = (
    #     world_length / 2 - margin_world_edge_source - d_source_det + 0.01 * m
    # )
    # create_phsp(
    #     sim,
    #     "phsp_front_det",
    #     "world",
    #     z_phsp_front_msc,
    #     m,
    #     ["KineticEnergy", "PrePositionLocal"],
    # )

    # Diamond detector definition
    # lateral shift of det to define which strip is facing the mb
    lat_shift = (232.5e-6 / 2 + 232.5e-6 * nb_offset_strips) * direction_offset * m

    # 136 strips * 0.1725 mm + 135 interstrips of 0.06 mm + 2 sides of diamond of 0.06mm
    # = 31.68 mm of detector width
    create_box_vol(
        sim,
        "diamond_detector",
        "world",
        [31.68e-3 * m, 3.32e-3 * m, 150e-6 * m],
        [lat_shift, 0, world_length / 2 - margin_world_edge_source - d_source_det],
        "diamant_det",
        [1, 0, 0, 1],
        cut_detector_photon,
        1 * m,
    )

    # volume containing all the strips for dose actor retrievements
    # 136 strips of 172.5 um and 135 interstrips of 60 um = 31.56 mm
    create_box_vol(
        sim,
        "dose_actor_vol",
        "diamond_detector",
        [31.56e-3 * m, 3.2e-3 * m, 150e-6 * m],
        [0, 0, 0],
        "diamant_det",
        [0, 0, 1, 1],
        cut_detector_photon,
        1 * m,
    )

    # force collision
    # force_coll = sim.add_actor("ForceCollisionActor", "force_coll")
    # force_coll.attached_to = "dose_actor_vol"

    # dose actor
    dose = sim.add_actor("DoseActor", "dose_detector")
    dose.attached_to = "dose_actor_vol"
    dose.output_filename = os.path.join(
        os.path.join(os.path.dirname(__file__), "output"),
        f"{slab_d}mm_{slab_material}_detector.mhd",
    )
    dose.dose.active = False
    dose.edep_uncertainty.active = True
    # pixel width = 0.0075 mm => 23 pixels per strips and 8 pixels per interstrip space
    # 23 * 136 strips + 8 * 135 interstrips = 4208 pixels
    dose.size = [4208, 1, 1]
    dose.spacing = [0.0075e-3 * m, 3.2e-3 * m, 150e-6 * m]

    # phys
    sim.physics_manager.physics_list_name = physics_list

    # stat actor
    s = sim.add_actor("SimulationStatisticsActor", "stats")
    s.track_types_flag = True
    s.output_filename = os.path.join(
        os.path.join(os.path.dirname(__file__), "output"), "stats.txt"
    )

    # ---------------------------------------------------------------------
    # start simulation
    if sleep_time:
        rd_time = np.random.randint(0, 120, 1)[0]
        time.sleep(rd_time)
    sim.run(start_new_process=False)

    stats = sim.get_actor("stats")
    print(stats)


@click.command()
@click.option("--n_part", default=1, help="Number of particle to simulate")
def main(n_part):
    run_simulation(n_part)


if __name__ == "__main__":
    main()
