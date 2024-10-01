import os
import sys
import opengate as gate
import numpy as np
import pandas as pd
import itk


# units
m = gate.g4_units.m
cm = gate.g4_units.cm
mm = gate.g4_units.mm
um = gate.g4_units.um
nm = gate.g4_units.nm
keV = gate.g4_units.keV


def create_array_of_elements(n_elements, center_to_center, *, offset_to_center=0):
    """Function creating an array of elements separated by a center_to_center distance and centered on 0 by default
    on the x-axis"""
    max_dist = center_to_center * (n_elements - 1)
    return np.linspace(-max_dist / 2, max_dist / 2, n_elements) + offset_to_center


def create_stripped_detector(
    mother_volume, n_strips, center_to_center, *, offset_to_center=0
):
    """Function creating a diamond detector with given number of strips separated by a center_to_center dist"""

    # gets positions of each strips center
    center_positions = create_array_of_elements(
        n_strips, center_to_center, offset_to_center=offset_to_center
    )

    for strip_number, pos in enumerate(center_positions):
        strip = sim.add_volume("BoxVolume", f"strip {strip_number}")
        strip.mother = mother_volume
        strip.size = [172.5e-6 * m, 3.2 * mm, 150e-6 * m]
        strip.material = "diamant_det"
        strip.color = [1, 0, 0, 1]
        strip.translation = [pos, 0, 0 * m]

        # dose actor
        dose = sim.add_actor("DoseActor", "dose_strip_" + str(strip_number))
        dose.mother = f"strip {strip_number}"
        dose.output = f"output/{rw3_thickness}cm_rw3_strip_{strip_number}.mhd"
        dose.square = False
        dose.dose = False
        dose.size = [1, 1, 1]

        # sets hit digitizers
        hits_digit = sim.add_actor(
            "DigitizerHitsCollectionActor", f"hits_{strip_number}"
        )
        hits_digit.mother = f"strip {strip_number}"
        hits_digit.output = f"output/{rw3_thickness}cm_rw3_strip_{strip_number}.root"
        hits_digit.attributes = ["TotalEnergyDeposit", "TrackID"]

        # sets production cut
        sim.physics_manager.set_production_cut(f"strip {strip_number}", "all", 0.1 * mm)


if __name__ == "__main__":
    # N_PARTICLES = 1_000_000_000
    # N_THREADS = 62
    N_PARTICLES = 10_000
    N_THREADS = 1
    rw3_thickness = 0

    # create the simulation
    sim = gate.Simulation()
    db_path = os.path.join(os.path.dirname(__file__), "data/gate_materials.db")
    sim.volume_manager.add_material_database(db_path)

    # main options
    sim.g4_verbose = False
    sim.visu = False
    sim.visu_type = "vrml"
    sim.check_volumes_overlap = False
    sim.number_of_threads = N_THREADS
    sim.random_seed = "auto"

    # world size
    sim.world.size = [1 * m, 1 * m, 8 * m]
    sim.world.material = "Air"

    # Microbeam definition
    source = sim.add_source("GenericSource", f"broad beam")
    source.particle = "gamma"
    source.n = N_PARTICLES / sim.number_of_threads
    source.position.type = "box"
    source.position.size = [50e-6 * m, 795e-6 * m, 1 * um]
    source.position.translation = [0, 0 * cm, -3.3 * m]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]
    source.energy.type = "mono"
    source.energy.mono = 123 * keV

    # RW3 slab definition
    if rw3_thickness != 0:
        wb = sim.add_volume("Box", "rw3_slab")
        wb.size = [10 * cm, 10 * cm, rw3_thickness * cm]
        wb.translation = [0, 0, -2 * m]
        wb.material = "RW3"
        wb.color = [0, 0, 1, 1]
        sim.physics_manager.set_production_cut("rw3_slab", "all", 1 * mm)

    # diamond detector definition
    diamond_detector = sim.add_volume("Box", "diamond_detector")
    # 8 diamonds of 4.013 mm width = 32.104 mm
    diamond_detector.size = [32.104 * mm, 3.32 * mm, 150e-6 * m]
    diamond_detector.material = "diamant_det"
    diamond_detector.color = [0, 1, 1, 1]

    # lateral shift of det
    # 136 strips (8 diamonds), if no lat shift mb centered between strip 68 and 69
    # to center mb on strip 60 (= strip 77 if 9 diamonds) : lat shift of 8 * pitch +
    # half strip widht and half interstrip distance
    # nb_offset_strips = 8
    # lat_shift = (nb_offset_strips * 232.5 + 172.5 / 2 + 0.06 / 2) * 10 ** (-6)
    lat_shift = 0
    diamond_detector.translation = [lat_shift, 0, 3.3 * m]
    sim.physics_manager.set_production_cut("diamond_detector", "all", 1 * mm)

    # strips detector definition
    # create_stripped_detector("diamond_detector", 136, 232.5e-6 * m)
    create_stripped_detector("diamond_detector", 17, 232.5e-6 * m)

    # phys
    sim.physics_manager.physics_list_name = "G4EmLivermorePolarizedPhysics"

    # stat actor
    s = sim.add_actor("SimulationStatisticsActor", "stats")
    s.track_types_flag = True
    s.output = os.path.join(
        os.path.join(os.path.dirname(__file__), "output"), "stats.txt"
    )

    # ---------------------------------------------------------------------
    # start simulation
    sim.run(start_new_process=False)

    stats = sim.output.get_actor("stats")
    print(stats)
