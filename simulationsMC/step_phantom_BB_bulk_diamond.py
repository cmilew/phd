import os
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


if __name__ == "__main__":

    N_PARTICLES = 100_000_000
    N_THREADS = 60
    water_thickness = 3

    # create the simulation
    sim = gate.Simulation()
    db_path = os.path.join(os.path.dirname(__file__), "additional_materials.db")
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
    sim.world.material = "G4_AIR"

    # Broad beam definition
    source = sim.add_source("GenericSource", f"broad beam")
    source.particle = "gamma"
    source.n = N_PARTICLES / sim.number_of_threads
    source.position.type = "box"
    source.position.size = [3 * cm, 3 * cm, 1 * nm]
    source.position.translation = [0, 0 * cm, -1 * m]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]
    source.energy.type = "mono"
    source.energy.mono = 123 * keV

    # Water slab definition
    wb = sim.add_volume("Box", "water_slab")
    wb.size = [10 * cm, 10 * cm, water_thickness * cm]
    wb.translation = [0, 0, 0 * m]
    wb.material = "G4_WATER"
    wb.color = [0, 0, 1, 1]  # blue
    sim.physics_manager.set_production_cut("water_slab", "all", 1 * mm)

    # Bulk diamond detector definition
    diamond_detector = sim.add_volume("Box", "diamond_detector")
    diamond_detector.size = [5 * cm, 5 * cm, 1 * cm]
    diamond_detector.material = "G4_C"
    diamond_detector.color = [0, 0, 1, 1]  
    diamond_detector.translation = [0, 0 * cm, 1 * m]
    sim.physics_manager.set_production_cut("diamond_detector", "all", 1 * mm)

    dose = sim.add_actor("DoseActor", "dose")
    dose.mother = "diamond_detector"
    dose.output = f"diamond_detector_{water_thickness}cm_water.mhd"
    # dose.spacing = [172.5e-6 * m, 3 * mm, 150e-6 * m]
    # dose.translation = [0, 0 * cm, 3.3 * m]
    # dose.uncertainty = True
    # dose.hit_type = "random"
    # dose.dose = True
    # dose.dose_calc_on_the_fly = (
    #     False  # calc dose as edep/mass after end of simulation
    # )
    dose.size = [1, 1, 1]

    # hits digitizer
    hits_digit = sim.add_actor('DigitizerHitsCollectionActor', "hits_diamond_detector")
    hits_digit.mother = "diamond_detector"
    path_output_file = os.path.join(os.path.dirname(__file__), f'diamond_detector{water_thickness}cm_water_hits.root')
    hits_digit.output = path_output_file
    hits_digit.attributes = ['TotalEnergyDeposit', 'TrackID']


    # phys
    sim.physics_manager.physics_list_name = "G4EmLivermorePolarizedPhysics"


    # stat actor
    s = sim.add_actor("SimulationStatisticsActor", "stats")
    s.track_types_flag = True
    s.output = os.path.join(os.path.join(os.path.dirname(__file__), "output") , "stats2.txt")

    # ---------------------------------------------------------------------
    # start simulation
    sim.run(start_new_process=True)

    stats = sim.output.get_actor("stats")
    print(stats)

        # read output
    d_post_path = sim.output.get_actor("dose").user_info.output
    img_mhd_out = itk.imread(d_post_path)
