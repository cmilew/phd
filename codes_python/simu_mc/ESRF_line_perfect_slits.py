import os
import sys
import opengate as gate
import numpy as np
import time
from opengate.geometry import volumes
from box import Box
from opengate.geometry.volumes import subtract_volumes, unite_volumes
import math
from decimal import Decimal


# units
m = gate.g4_units.m
cm = gate.g4_units.cm
mm = gate.g4_units.mm
um = gate.g4_units.um
nm = gate.g4_units.nm
keV = gate.g4_units.keV
MeV = gate.g4_units.MeV
deg = gate.g4_units.deg


def create_box_vol(
    box_name, box_mother, box_size, box_translation, box_material, box_color
):
    """Function creating a box volume and setting its properties"""

    box = volumes.BoxVolume(name=box_name)
    box.mother = box_mother
    box.size = box_size
    box.translation = box_translation
    box.material = box_material
    box.color = box_color

    return box


def create_collimator(
    mother_volume,
    col_box_name,
    col_name,
    slit_name,
    material,
    translation,
    box_dim,
    slit_dim,
    cut_col,
    cut_slit,
):
    """Function creating a collimator by a subtraction of a slit to a box, adds
    collimator and its slit to simulation and set their cuts"""

    # Create box
    col_box = create_box_vol(
        col_box_name, mother_volume, box_dim, translation, material, [1, 0, 0, 1]
    )

    # Create slit
    slit = create_box_vol(
        slit_name, col_box_name, slit_dim, [0, 0, 0], "G4_AIR", [0, 0, 0, 0]
    )

    # Create collimator = box - slit
    col = subtract_volumes(col_box, slit, new_name=col_name)

    # Add volumes to simulation
    # Mandatory to col box to add slit to the simulation (and slit needed to adjust cut)
    sim.add_volume(col_box)
    sim.add_volume(col)
    sim.add_volume(slit)
    sim.physics_manager.set_production_cut(col_name, "all", cut_col)
    sim.physics_manager.set_production_cut(slit_name, "all", cut_slit)


def create_array_of_elements(n_elements, center_to_center, *, offset_to_center=0):
    """Function creating an array of elements separated by a center_to_center distance
    and centered on 0 by default on the x-axis"""

    max_dist = center_to_center * (n_elements - 1)
    return np.linspace(-max_dist / 2, max_dist / 2, n_elements) + offset_to_center


def create_msc(
    d_source_msc,
    n_slits,
    center_to_center,
    cut_msc,
    cut_slits_msc,
    *,
    offset_to_center=0,
):
    """Function creating the MSC as a subtraction of the slits to a box, unite every
    slit of MSC, adds MSC and united slits to simulation and set their cuts"""

    # Create MSC box
    msc_box = create_box_vol(
        "msc_box",
        "world",
        [8 * cm, 1 * cm, 8 * mm],
        [0, 0, d_source_msc * m],
        "msc_material",
        [1, 0, 0, 1],
    )

    # Gets positions of each slit center
    center_positions = create_array_of_elements(
        n_slits, center_to_center, offset_to_center=offset_to_center
    )
    msc = msc_box

    # Initialize msc_slits volume with the middle slit (= slit 63 as they are 125 slits)
    msc_slits = create_box_vol(
        "msc_slit_0",
        "msc_box",
        [50 * um, 3 * mm, 8.1 * mm],
        [center_positions[len(center_positions) // 2], 0, 0 * m],
        "G4_AIR",
        [0, 0, 0, 0],
    )

    # Subtract slits to MSC box and unite them
    for slit_num, pos in enumerate(center_positions):
        slit = volumes.BoxVolume(name=f"slit {slit_num}")
        slit.mother = msc_box
        slit.size = [50 * um, 3 * mm, 8.1 * mm]
        slit.material = "G4_AIR"
        msc = subtract_volumes(msc, slit, translation=[pos, 0, 0 * m], new_name="msc")
        msc_slits = unite_volumes(
            msc_slits, slit, translation=[pos, 0, 0 * m], new_name="msc_slits"
        )

    sim.add_volume(msc_box)
    sim.add_volume(msc)
    sim.add_volume(msc_slits)
    sim.physics_manager.set_production_cut("msc", "all", cut_msc)
    sim.physics_manager.set_production_cut("msc", "all", cut_slits_msc)


if __name__ == "__main__":
    # Simulation parameters
    N_PARTICLES = 10
    N_THREADS = 1
    sleep_time = False  # only for parallel simulation on CC
    physics_list = "G4EmLivermorePolarizedPhysics"
    phsp = True
    visu = True
    kill_act = True

    # Maximale distance that need to be covered by the beam (half diagonal of detector
    #  if full detector needs to be irradiated) => will define divergence of source (theta)
    # for having at least 800 um beam height (a little more than 795 um) at vertical slit dist
    max_coverage_dist = 400e-6
    # Distance between source and where the max coverage dist is calculated
    d_target_source = 40.2
    beam_div = np.arctan(max_coverage_dist / d_target_source) * (180 / math.pi)
    print(beam_div)
    beam_height = 795 * um

    prim_col_cut = 1 * m
    prim_slit_cut = 1 * mm
    sec_col_cut = 1 * mm
    sec_slit_cut = 10 * um
    msc_cut = 10 * um  # short cut to be precise in MSC leaves
    msc_slits_cut = 10 * um

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
    # theta chosen to have a beam width of one microbeam width at MSC
    source.direction.theta = [0 * deg, beam_div * deg]
    source.direction.phi = [0 * deg, 360 * deg]
    source.energy.type = "spectrum_lines"
    spectrum = np.loadtxt(
        os.path.join(os.path.dirname(__file__), "data/ESRF_clinic_spec.txt"),
        skiprows=1,
        delimiter="\t",
    )
    source.energy.spectrum_energy = spectrum[:, 0] * MeV
    source.energy.spectrum_weight = spectrum[:, 1]

    # Primary collimator, slit of 43 um width to get beam width of 60 um at MSC
    # (a little more than one microbeam width) and 583 um height to get 800 um at sec col
    # (a little more than 795 um)
    create_collimator(
        "vacuum_sec",
        "pc_box",
        "primary_collimator",
        "primary_slit",
        "Copper",
        [0, 0, -10.75 * m],
        [5 * cm, 5 * cm, 1 * cm],
        [43 * um, 583 * um, 1.1 * cm],
        prim_col_cut,
        prim_slit_cut,
    )

    # Secondary collimator, slit of 59 um width to get beam width of 60 um at MSC
    # beam height defined by user
    create_collimator(
        "world",
        "sc_box",
        "secondary_collimator",
        "secondary_slit",
        "msc_material",
        [0, 0, -19.2 * m],
        [8 * cm, 1 * cm, 8 * mm],
        [59 * um, beam_height, 8.1 * mm],
        sec_col_cut,
        sec_slit_cut,
    )

    # Creates MSC
    msc = create_msc(-19.7, 125, 400 * um, msc_cut, msc_slits_cut, offset_to_center=0)

    # kill actors
    if kill_act:
        kill_actor_pc = sim.add_actor("KillActor", "kill_actor_pc")
        kill_actor_pc.mother = "primary_collimator"
        kill_actor_sc = sim.add_actor("KillActor", "kill_actor_sc")
        kill_actor_sc.mother = "secondary_collimator"
        kill_actor_mlc = sim.add_actor("KillActor", "kill_actor_msc")
        kill_actor_mlc.mother = "msc"

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
            f"output/phsp_esrf_line_{Decimal(N_PARTICLES):.3E}_events.root",
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
    if sleep_time:
        rd_time = np.random.randint(0, 120, 1)[0]
        time.sleep(rd_time)
    sim.run(start_new_process=False)

    stats = sim.output.get_actor("stats")
    print(stats)
