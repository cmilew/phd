import os
import sys
import opengate as gate
import numpy as np
import pandas as pd
import itk
import uproot
import click
import time


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


def create_array_of_elements(n_elements, center_to_center, *, offset_to_center=0):
    """Function creating an array of elements separated by a center_to_center distance and centered on 0 by default
    on the x-axis"""
    max_dist = center_to_center * (n_elements - 1)
    return np.linspace(-max_dist / 2, max_dist / 2, n_elements) + offset_to_center


def create_dose_actor_strips(
    sim,
    slab_thickness,
    m,
    mother_volume,
    n_strips,
    center_to_center,
    *,
    offset_to_center=0,
):
    """Function creating a dose actor for every strips each of them separated by a center_to_center dist"""

    # gets positions of each strips center
    center_positions = create_array_of_elements(
        n_strips, center_to_center, offset_to_center=offset_to_center
    )

    for strip_number, pos in enumerate(center_positions):
        strip = sim.add_volume("BoxVolume", f"strip {strip_number}")
        strip.mother = mother_volume
        strip.size = [172.5e-6 * m, 3.2e-3 * m, 150e-6 * m]
        strip.material = "diamant_det"
        strip.color = [1, 0, 0, 1]
        strip.translation = [pos, 0, 0 * m]

        # dose actor
        dose = sim.add_actor("DoseActor", "dose_strip_" + str(strip_number))
        dose.mother = f"strip {strip_number}"
        dose.output = f"output/{slab_thickness}cm_rw3_strip_{strip_number}.mhd"
        dose.square = False
        dose.dose = False
        dose.size = [1, 1, 1]

        # sets hit digitizers
        hits_digit = sim.add_actor(
            "DigitizerHitsCollectionActor", f"hits_{strip_number}"
        )
        hits_digit.mother = f"strip {strip_number}"
        hits_digit.output = f"output/{slab_thickness}cm_rw3_strip_{strip_number}.root"
        hits_digit.attributes = ["TotalEnergyDeposit", "TrackID"]

        # sets production cut
        sim.physics_manager.set_production_cut(
            f"strip {strip_number}", "all", 0.1e-3 * m
        )


def run_simulation(n_part):
    # units
    m = gate.g4_units.m

    # simulation parameters
    N_PARTICLES = 100_000
    N_THREADS = 1
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
    phsp_path = os.path.join(os.path.dirname(__file__), f"data/{phsp_filename}")
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
    # strips detector definition
    create_stripped_detector(sim, slab_d, m, "diamond_detector", 136, 232.5e-6 * m)

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
