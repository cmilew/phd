import os
import sys
import opengate as gate
import numpy as np
import time
from opengate.geometry import volumes
from opengate.geometry.volumes import intersect_volumes
from decimal import Decimal

# units
m = gate.g4_units.m
cm = gate.g4_units.cm
mm = gate.g4_units.mm
um = gate.g4_units.um
nm = gate.g4_units.nm
keV = gate.g4_units.keV
deg = gate.g4_units.deg


def create_array_of_elements(n_elements, center_to_center, *, offset_to_center=0):
    """Function creating an array of elements separated by a center_to_center distance and centered on 0 by default
    on the x-axis"""
    max_dist = center_to_center * (n_elements - 1)
    return np.linspace(-max_dist / 2, max_dist / 2, n_elements) + offset_to_center


def create_mlc_slits(mlc_vol, n_slits, center_to_center, *, offset_to_center=0):
    """Function creating a the slits in MLC volume"""

    # gets positions of each strips center
    center_positions = create_array_of_elements(
        n_slits, center_to_center, offset_to_center=offset_to_center
    )

    for slit_num, pos in enumerate(center_positions):
        slit = sim.add_volume("BoxVolume", f"slit {slit_num}")
        slit.mother = mlc_vol
        slit.size = [50 * um, 3 * mm, 8 * mm]
        slit.material = "G4_AIR"
        slit.translation = [pos, 0, 0 * m]
        sim.physics_manager.set_production_cut(f"slit {slit_num}", "all", 10 * um)


if __name__ == "__main__":
    N_PARTICLES = 1_000_000
    N_THREADS = 18
    energy = 123
    physics_list = "G4EmLivermorePolarizedPhysics"
    phsp = True
    visu = False

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
    sim.world.size = [1 * m, 1 * m, 42 * m]
    sim.world.material = "Air"
    sim.physics_manager.set_production_cut("world", "all", 1 * m)

    # vacuum section
    vac_sec = sim.add_volume("Box", "vacuum_sec")
    vac_sec.size = [1 * m, 1 * m, 37.1 * m]
    vac_sec.translation = [0, 0, 2.45 * m]
    vac_sec.material = "Vacuum"
    vac_sec.color = [0, 0.5, 0, 0.5]
    sim.physics_manager.set_production_cut("vacuum_sec", "all", 1 * m)

    # Point source defined at wiggler
    source = sim.add_source("GenericSource", "point_source")
    source.particle = "gamma"
    source.mother = "world"
    source.n = N_PARTICLES / sim.number_of_threads
    source.position.type = "point"
    source.position.translation = [0, 0, 21 * m]
    source.direction.type = "iso"
    source.direction.theta = [0 * deg, 3.52e-5 * deg]
    source.direction.phi = [0 * deg, 360 * deg]
    source.energy.type = "mono"
    source.energy.mono = energy * keV

    # Primary collimator
    prim_col = sim.add_volume("Box", "primary_collimator")
    prim_col.mother = "vacuum_sec"
    prim_col.size = [5 * cm, 5 * cm, 1 * cm]
    prim_col.translation = [0, 0, -10.75 * m]
    prim_col.material = "Copper"
    prim_col.color = [1, 0, 0, 1]
    sim.physics_manager.set_production_cut("primary_collimator", "all", 1 * m)

    # Primary slit in primary collimator, define beam width and height
    prim_slit = sim.add_volume("Box", "prim_slit")
    # one microbeam width aperture
    prim_slit.size = [36 * um, 795 * um, 1.1 * cm]
    prim_slit.material = "G4_AIR"
    prim_slit.mother = "primary_collimator"
    prim_slit.translation = [0, 0, 0]
    sim.physics_manager.set_production_cut("prim_slit", "all", 10 * um)

    # MLC
    mlc = sim.add_volume("Box", "mlc")
    mlc.mother = "world"
    mlc.size = [8 * cm, 1 * cm, 8 * mm]
    mlc.translation = [0, 0, -19.7 * m]
    mlc.material = "mlc_material"
    create_mlc_slits("mlc", 125, 400 * um, offset_to_center=0)
    sim.physics_manager.set_production_cut("mlc", "all", 1 * mm)

    # # kill actors
    # kill_actor_pc = sim.add_actor("KillActor", "kill_actor_pc")
    # kill_actor_pc.attached_to = "primary_collimator"
    # kill_actor_mlc = sim.add_actor("KillActor", "kill_actor_mlc")
    # kill_actor_mlc.attached_to = "mlc"

    # virtual plane for phase space
    if phsp:
        plane = sim.add_volume("Box", "phsp_plane")
        plane.mother = sim.world
        plane.material = "G4_AIR"
        plane.size = [20 * cm, 10 * cm, 1 * um]
        plane.translation = [0, 0, -20 * m]
        plane.color = [1, 0, 0, 1]

        # phaseSpace Actor
        phsp_actor = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
        phsp_actor.mother = plane.name
        phsp_actor.attributes = [
            "KineticEnergy",
            "Weight",
            "PrePositionLocal",
            "PreDirectionLocal",
        ]
        phsp_actor.output = os.path.join(
            os.path.dirname(__file__),
            f"output/phsp_esrf_{energy}kev_{Decimal(N_PARTICLES):.3E}_events.root",
        )
        phsp_actor.debug = False

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
    # rd_time = np.random.randint(0, 120, 1)[0]
    # time.sleep(rd_time)
    sim.run(start_new_process=False)

    stats = sim.output.get_actor("stats")
    print(stats)
