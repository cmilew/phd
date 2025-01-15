import os
import opengate as gate


if __name__ == "__main__":
    # units
    m = gate.g4_units.m
    keV = gate.g4_units.keV
    deg = gate.g4_units.deg

    # simulation parameters
    N_PARTICLES = 10
    N_THREADS = 1
    visu = False
    physics_list = "G4EmLivermorePolarizedPhysics"

    world_length = 6 * m

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

    # source parameters
    source = sim.add_source("GenericSource", "point_source")
    source.particle = "gamma"
    source.attached_to = "world"
    source.n = N_PARTICLES / sim.number_of_threads
    source.position.type = "box"
    source.position.size = [0.05 * m, 0.05 * m, 0.01 * m]
    source.position.translation = [0, 0, 0]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]
    source.energy.type = "mono"
    source.energy.mono = 106 * keV

    # phys
    sim.physics_manager.physics_list_name = physics_list

    s = f"/process/em/UseGeneralProcess false"
    sim.g4_commands_before_init.append(s)

    # ---------------------------------------------------------------------
    # start simulation
    sim.run()
