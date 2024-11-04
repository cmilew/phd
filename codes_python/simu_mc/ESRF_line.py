import os
import opengate as gate
import numpy as np
import time
from opengate.geometry import volumes
from opengate.geometry.volumes import subtract_volumes, unite_volumes
import math
from decimal import Decimal
import click


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
    simulation,
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
    # Mandatory to add col box to add slit to the simulation (and slit needed to adjust cut)
    simulation.add_volume(col_box)
    simulation.add_volume(col)
    simulation.add_volume(slit)
    simulation.physics_manager.set_production_cut(col_name, "all", cut_col)
    simulation.physics_manager.set_production_cut(slit_name, "all", cut_slit)


def create_array_of_elements(n_elements, center_to_center, *, offset_to_center=0):
    """Function creating an array of elements separated by a center_to_center distance
    and centered on 0 by default on the x-axis"""

    max_dist = center_to_center * (n_elements - 1)
    return np.linspace(-max_dist / 2, max_dist / 2, n_elements) + offset_to_center


def create_msc(
    simulation,
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

    m = gate.g4_units.m
    mm = gate.g4_units.mm

    # Create MSC box
    msc_box = create_box_vol(
        "msc_box",
        "world",
        [80 * mm, 10 * mm, 8 * mm],
        [0, 0, d_source_msc],
        "msc_material",
        [1, 0, 0, 1],
    )

    # Gets positions of each slit center
    center_positions = create_array_of_elements(
        n_slits, center_to_center, offset_to_center=offset_to_center
    )
    msc = msc_box

    # Initialize msc_slits volume (for unite boolean) with the middle slit
    # (= slit 63 as they are 125 slits)
    msc_slits = create_box_vol(
        "msc_slit_0",
        "msc_box",
        [50e-3 * mm, 3 * mm, 8.1 * mm],
        [center_positions[len(center_positions) // 2], 0, 0 * m],
        "G4_AIR",
        [0, 0, 0, 0],
    )

    # Subtract slits to MSC box and unite them
    for slit_num, pos in enumerate(center_positions):
        slit = volumes.BoxVolume(name=f"slit {slit_num}")
        slit.mother = msc_box
        slit.size = [50e-3 * mm, 3 * mm, 8.1 * mm]
        slit.material = "G4_AIR"
        msc = subtract_volumes(msc, slit, translation=[pos, 0, 0 * m], new_name="msc")
        msc_slits = unite_volumes(
            msc_slits, slit, translation=[pos, 0, 0 * m], new_name="msc_slits"
        )

    simulation.add_volume(msc_box)
    simulation.add_volume(msc)
    simulation.add_volume(msc_slits)
    simulation.physics_manager.set_production_cut("msc", "all", cut_msc)
    simulation.physics_manager.set_production_cut("msc", "all", cut_slits_msc)


def run_simulation(n_part):
    # units
    m = gate.g4_units.m
    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    um = gate.g4_units.um
    eV = gate.g4_units.eV
    keV = gate.g4_units.keV
    deg = gate.g4_units.deg

    # Simulation parameters
    N_PARTICLES = n_part
    N_THREADS = 18
    sleep_time = False  # only for parallel simulation on CC
    physics_list = "G4EmLivermorePolarizedPhysics"
    phsp = True
    visu = False
    kill_act = False

    # Beam dimensions defined at sec col
    # beam_width = 35e-3 * m  # beam width for step phantom setup
    beam_width = 60e-6 * m  # 1 mb
    beam_height = 795e-6 * m

    # Dist between source and collimators
    d_source_prim_col = 29.3 * m
    d_source_sec_col = 40.2 * m

    # Beam divergence calc so that beam width = half diagonal of sec col slit + margin
    # (to cover the whole slit of sec col)
    margin = 10e-6 * m
    hald_diag_sec_col = (1 / 2) * np.sqrt(beam_width**2 + beam_height**2)
    beam_div = (
        np.arctan((hald_diag_sec_col + margin) / d_source_sec_col)
        * (180 / math.pi)
        * deg
    )

    # Cuts
    prim_col_cut = 1 * m
    prim_slit_cut = 1 * mm
    sec_col_cut = 1 * mm
    sec_slit_cut = 1 * mm
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
    source.direction.theta = [0 * deg, beam_div]
    source.direction.phi = [0 * deg, 360 * deg]
    source.energy.type = "mono"
    source.energy.mono = 123 * keV
    # source.energy.type = "spectrum_lines"
    # spectrum = np.loadtxt(
    #     os.path.join(os.path.dirname(__file__), "data/ESRF_clinic_spec_maxi_bins.txt"),
    #     skiprows=1,
    #     delimiter="\t",
    # )
    # source.energy.spectrum_energy = spectrum[:, 0] * eV
    # source.energy.spectrum_weight = spectrum[:, 1]

    # Primary collimator
    # prim col width chosen to cover whole sec col width + margin on both sides
    prim_col_width = (beam_width + 2 * margin) * (d_source_prim_col / d_source_sec_col)
    prim_col_height = (beam_height + 2 * margin) * (
        d_source_prim_col / d_source_sec_col
    )
    create_collimator(
        sim,
        "vacuum_sec",
        "pc_box",
        "primary_collimator",
        "primary_slit",
        "Copper",
        [0, 0, -10.75 * m],
        [5 * cm, 5 * cm, 1 * cm],
        [prim_col_width, prim_col_height, 1.1 * cm],
        prim_col_cut,
        prim_slit_cut,
    )

    # Secondary collimator
    create_collimator(
        sim,
        "world",
        "sc_box",
        "secondary_collimator",
        "secondary_slit",
        "msc_material",
        [0, 0, -19.2 * m],
        [8 * cm, 1 * cm, 8 * mm],
        [beam_width, beam_height, 8.1 * mm],
        sec_col_cut,
        sec_slit_cut,
    )

    # Creates MSC
    create_msc(
        sim, -19.7 * m, 125, 400 * um, msc_cut, msc_slits_cut, offset_to_center=0
    )

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
            "ParticleName",
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


@click.command()
@click.option("--n_part", default=1, help="Number of particle to simulate")
def main(n_part):
    run_simulation(n_part)


if __name__ == "__main__":
    main()
