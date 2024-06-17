#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import opengate as gate
import numpy as np
import pandas as pd
import uproot


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


def create_microbeam_array(n_microbeams, center_to_center, *, offset_to_center=0):
    """Function creating an array of n microbeams, each center separated by a center_to_center dist to another"""
    # gets positions of each microbeam center
    center_positions = create_array_of_elements(n_microbeams, center_to_center, offset_to_center=offset_to_center)

    # retrieve microbeams spectrum
    spectrum_path = os.path.join(os.path.dirname(__file__), "ESRF_spectrum.xlsx")
    df = pd.read_excel(spectrum_path, skiprows=6)
    spectrum_energy = df['Energy [eV]'].to_numpy()
    spectrum_weight = df['Flux_1'].to_numpy()
    spectrum_weight = spectrum_weight / np.sum(spectrum_weight)

    for beam_number, pos in enumerate(center_positions):
        source = sim.add_source("GenericSource", f"microbeam {beam_number}")
        source.particle = "gamma"
        source.n = 100_000 / sim.number_of_threads
        source.position.type = "box"
        source.position.size = [50e-6 * m, 795e-6 * m, 1 * nm]
        source.position.translation = [pos, 0 * cm, -3.3 * m]
        source.direction.type = "momentum"
        source.direction.momentum = [0, 0, 1]
        source.energy.type = "spectrum_lines"
        source.energy.spectrum_weight = spectrum_weight
        source.energy.spectrum_energy = spectrum_energy * keV


def create_stripped_detector(n_strips, center_to_center, *, offset_to_center=0):
    """Function creating a diamond detector with given number of strips separated by a center_to_center dist"""
    # gets positions of each strips center
    center_positions = create_array_of_elements(n_strips, center_to_center, offset_to_center=offset_to_center)

    for strip_number, pos in enumerate(center_positions):
        strip = sim.add_volume("BoxVolume", f"strip {strip_number}")
        strip.size = [172.5e-4 * m, 3 * mm, 150e-6 * m]
        strip.material = "G4_C"
        strip.color = [220, 38, 127, 1]
        strip.translation = [pos, 0 * cm, 3.3 * m]

        # sets hit digitizers
        hits_digit = sim.add_actor('DigitizerHitsCollectionActor', f"hits {strip_number}")
        hits_digit.mother = f"strip_{strip_number}"
        path_output_file = os.path.join(os.path.dirname(__file__), f'{strip_number}_hits.root')
        hits_digit.output = path_output_file
        hits_digit.attributes = ['TotalEnergyDeposit', 'TrackID']

        # sets production cut
        sim.physics_manager.set_production_cut(f"strip {strip_number}", "all", 0.1 * mm)



if __name__ == "__main__":

    # create the simulation
    sim = gate.Simulation()
    db_path = os.path.join(os.path.dirname(__file__), "additional_materials.db")
    sim.volume_manager.add_material_database(db_path)

    # main options
    sim.g4_verbose = False
    sim.visu = False
    sim.visu_type = "vrml"
    sim.check_volumes_overlap = False
    sim.number_of_threads = 1
    sim.random_seed = "auto"

    # world size
    sim.world.size = [1 * m, 1 * m, 8 * m]
    sim.world.material = "G4_AIR"

    # RW3 slab definition
    wb = sim.add_volume("Box", "RW3_slab")
    wb.size = [16 * cm, 2 * cm, 1 * cm]
    wb.translation = [0, 0, -2 * m]
    wb.material = "RW3"
    wb.color = [0, 0, 1, 1]  # blue
    wb.color = [1, 0.7, 0.7, 0.8]

    # PBC part of the detector definition
    # PBC_detector = sim.add_volume("Box", "PVC_detector")
    # PBC_detector.size = [20 * cm, 10 * cm, 1 * cm]
    # PBC_detector.material = "PBC"
    # PBC_detector.color = [0.5, 0.5, 0.5, 1]
    # PBC_detector.translation = [0, 0 * cm, 3.3 * m]

    # microbeam array definition (shifted of 200µm to align 25th microbeam with central strip of detector)
    create_microbeam_array(50, 400e-6 * m, offset_to_center=200e-6 * m)

    # detector definition
    # create_stripped_detector(5, 232.5e-4 * m)
    central_strip = sim.add_volume("Box", "central_strip")
    central_strip.size = [172.5e-6 * m, 3 * mm, 150e-6 * m]
    central_strip.material = "G4_C"
    central_strip.color = [0, 0, 1, 1]  
    central_strip.translation = [0, 0 * cm, 3.3 * m]
    sim.physics_manager.set_production_cut("central_strip", "all", 2.5 * um)
    hits_digit = sim.add_actor('DigitizerHitsCollectionActor', "hits_central_strip")
    hits_digit.mother = "central_strip"
    path_output_file = os.path.join(os.path.dirname(__file__), 'central_strip_hits.root')
    hits_digit.output = path_output_file
    hits_digit.attributes = ['TotalEnergyDeposit', 'TrackID']


    # phys
    sim.physics_manager.physics_list_name = "G4EmLivermorePolarizedPhysics"
    sim.physics_manager.set_production_cut("world", "all", 10 * m)
    sim.physics_manager.set_production_cut("RW3_slab", "all", 1 * mm)

    # stat actor
    s = sim.add_actor("SimulationStatisticsActor", "stats")
    s.track_types_flag = True
    s.output = os.path.join(os.path.join(os.path.dirname(__file__), "output") , "stats2.txt")

    # ---------------------------------------------------------------------
    # start simulation
    sim.run(start_new_process=True)

    stats = sim.output.get_actor("stats")
    print(stats)
