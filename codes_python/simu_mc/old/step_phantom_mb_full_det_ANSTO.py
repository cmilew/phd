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
        dose.output = f"output/{slab_depth}cm_rw3_strip_{strip_number}.mhd"
        dose.square = False
        dose.dose = False
        dose.size = [1, 1, 1]

        # sets hit digitizers
        hits_digit = sim.add_actor(
            "DigitizerHitsCollectionActor", f"hits_{strip_number}"
        )
        hits_digit.mother = f"strip {strip_number}"
        hits_digit.output = f"output/{slab_depth}cm_rw3_strip_{strip_number}.root"
        hits_digit.attributes = ["TotalEnergyDeposit", "TrackID"]

        # sets production cut
        sim.physics_manager.set_production_cut(f"strip {strip_number}", "all", 0.1 * mm)


if __name__ == "__main__":
    N_PARTICLES = 1
    N_THREADS = 1
    nb_offset_strips = 0
    # negative = towards 136, position = towards 1
    direction_offset = 1
    energy = 106
    slab_material = "SolidHE"
    slab_depth = 6
    slab_l = 20 * cm
    slab_L = 20 * cm
    physics_list = "G4EmStandardPhysics_option3"
    # physics_list = "G4EmLivermorePolarizedPhysics"

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
    sim.world.size = [1 * m, 1 * m, 10 * m]
    sim.world.material = "Air"
    # sim.physics_manager.set_production_cut("world", "all", 10 * m)

    # Microbeam definition
    source = sim.add_source("GenericSource", "broad beam")
    source.particle = "gamma"
    source.n = N_PARTICLES / sim.number_of_threads
    source.position.type = "box"
    source.position.size = [50 * um, 1.056 * mm, 1 * um]
    source.position.translation = [0, 0 * cm, -50 * cm]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]
    source.energy.type = "mono"
    source.energy.mono = energy * keV

    # slab definition
    if slab_depth != 0:
        slab = sim.add_volume("Box", "slab")
        slab.size = [slab_L * cm, slab_l * cm, slab_depth * cm]
        slab.translation = [0, 0, 0 * m]
        slab.material = slab_material
        slab.color = [0, 0, 1, 1]
        sim.physics_manager.set_production_cut("slab", "all", 1 * mm)

    # virtual plane for phase space
    plane = sim.add_volume("Box", "phsp_plane")
    plane.mother = sim.world
    plane.material = "G4_AIR"
    plane.size = [10 * cm, 5 * cm, 1 * um]
    plane.translation = [0, 0, 4 * m]
    plane.color = [1, 0, 0, 1]

    # phaseSpace Actor
    phsp_actor = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
    phsp_actor.mother = plane.name
    phsp_actor.attributes = ["KineticEnergy", "EventPosition"]
    phsp_actor.output = os.path.join(
        os.path.dirname(__file__),
        f"output/phsp_ansto_{energy}kev_{N_PARTICLES}_events.root",
    )
    phsp_actor.debug = False

    # diamond detector definition
    diamond_detector = sim.add_volume("Box", "diamond_detector")
    diamond_detector.mother = sim.world
    # 8 diamonds * 4.013 mm width = 32.104 mm
    diamond_detector.size = [32.104 * mm, 3.32 * mm, 150e-6 * m]
    diamond_detector.material = "diamant_det"
    diamond_detector.color = [0, 1, 1, 1]

    # lateral shift of det to define which strip is facing the mb
    lat_shift = (232.5 / 2 + 232.5 * nb_offset_strips) * direction_offset
    diamond_detector.translation = [lat_shift * um, 0, 4.13 * m]
    sim.physics_manager.set_production_cut("diamond_detector", "all", 1 * um)

    # strips detector definition
    create_stripped_detector("diamond_detector", 136, 232.5e-6 * m)

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
    sim.run(start_new_process=False)

    stats = sim.output.get_actor("stats")
    print(stats)

    from opengate.managers import ActorManager

    print(isinstance(sim.actor_manager, ActorManager))
    for i in sim.actor_manager.user_info_actors.values():
        print(type(i))
