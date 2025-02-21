// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <string>
#include <vector>
#include "Framework/Variant.h"
#include "Framework/ConfigParamSpec.h"
#include "DataFormatsEMCAL/Digit.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "EMCALWorkflow/PublisherSpec.h"
#include "CommonUtils/ConfigurableParam.h"

using namespace o2::framework;
using namespace o2::emcal;

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"disable-mc", VariantType::Bool, false, {"Do not propagate MC labels"}},
                                       {"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings"}}};
  o2::raw::HBFUtilsInitializer::addConfigOption(options);
  workflowOptions.insert(workflowOptions.end(), options.begin(), options.end());
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  bool disableMC = cfgc.options().get<bool>("disable-mc");
  o2::conf::ConfigurableParam::updateFromString(cfgc.options().get<std::string>("configKeyValues"));

  WorkflowSpec specs;
  specs.emplace_back(o2::emcal::getPublisherSpec<std::vector<o2::emcal::Digit>>(PublisherConf{
                                                                                  "emcal-digit-reader",
                                                                                  "o2sim",
                                                                                  "emcdigits.root",
                                                                                  {"digitbranch", "EMCALDigit", "Digit branch"},
                                                                                  {"digittriggerbranch", "EMCALDigitTRGR", "Trigger record branch"},
                                                                                  {"mcbranch", "EMCALDigitMCTruth", "MC label branch"},
                                                                                  o2::framework::OutputSpec{"EMC", "DIGITS"},
                                                                                  o2::framework::OutputSpec{"EMC", "DIGITSTRGR"},
                                                                                  o2::framework::OutputSpec{"EMC", "DIGITSMCTR"}},
                                                                                !disableMC));

  // configure dpl timer to inject correct firstTFOrbit: start from the 1st orbit of TF containing 1st sampled orbit
  o2::raw::HBFUtilsInitializer hbfIni(cfgc, specs);

  return specs;
}