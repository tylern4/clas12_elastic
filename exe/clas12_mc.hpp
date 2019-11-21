/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "histogram.hpp"
#include "reaction.hpp"

size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram>& _hists, int thread_id);
size_t run_files(std::vector<std::string> inputs, const std::shared_ptr<Histogram>& hists, int thread_id);

size_t run_files(std::vector<std::string> inputs, const std::shared_ptr<Histogram>& hists, int thread_id) {
  // Called once for each thread
  // Make a new chain to process for this thread
  auto chain = std::make_shared<TChain>("clas12");
  // Add every file to the chain
  for (auto in : inputs) chain->Add(in.c_str());

  // Run the function over each thread
  return run(chain, hists, thread_id);
}

size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram>& _hists, int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();
  float beam_energy = NAN;
  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  auto data = std::make_shared<Branches12>(_chain, true);

  // Total number of events "Processed"
  size_t total = 0;
  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);
    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    if (data->mc_npart() < 1) continue;

    // If we pass electron cuts the event is processed
    total++;

    // Make a reaction class from the data given
    auto mc_event = std::make_shared<MCReaction>(data, beam_energy);
    for (int part = 1; part < data->mc_npart(); part++) {
      // Check particle ID's and fill the reaction class
      if (data->mc_pid(part) == PIP) {
        mc_event->SetPip(part);
      } else if (data->mc_pid(part) == PROTON) {
        mc_event->SetProton(part);
      } else if (data->mc_pid(part) == PIM) {
        mc_event->SetPim(part);
      } else {
        mc_event->SetOther(part);
      }
    }
    _hists->Fill_WvsQ2(mc_event);

    // Assume mc event is regular reconstructed event
    if (data->gpart() == 0) continue;
    bool elec = true;
    elec &= (data->charge(0) == NEGATIVE);
    elec &= (data->pid(0) == 11);
    if (!elec) continue;

    auto event = std::make_shared<Reaction>(data, beam_energy);
    auto dt = std::make_shared<Delta_T>(data);
    // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);
      _hists->Fill_MomVsBeta(data, part);
      _hists->Fill_deltat_pi(data, dt, part);
      _hists->Fill_deltat_prot(data, dt, part);

      // Check particle ID's and fill the reaction class
      if (abs(dt->dt_Pi()) < 0.5 && data->charge(part) == POSITIVE) {
        event->SetPip(part);
      } else if (abs(dt->dt_P()) < 0.5 && data->charge(part) == POSITIVE) {
        event->SetProton(part);
      } else if (abs(dt->dt_Pi()) < 0.5 && data->charge(part) == NEGATIVE) {
        event->SetPim(part);
      } else {
        event->SetOther(part);
      }
    }

    // Check the reaction class what kind of even it is and fill the appropriate histograms
    if (event->SingleP()) _hists->Fill_WvsQ2_singleP(event);
    if (event->NeutronPip()) _hists->Fill_WvsQ2_Npip(event);
  }

  // Return the total number of events
  return num_of_events;
}
#endif
