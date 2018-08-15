#include "CentrifugeConfig.h"
#include "StageConfig.h"
#include "CascadeConfig.h"
#include <iostream>


int main () {
 

  int secpertimestep = 2629846;


  double design_feed_assay = 0.0071;
  double design_product_assay = 0.035;
  double design_tails_assay = 0.003;
  
  double design_feed_flow = 1000;
  int max_centrifuges = 169;
  double precision = 1.e-15;

  mbmore::CentrifugeConfig centrifuge = mbmore::CentrifugeConfig();
  // Update Centrifuge paramter from the user input:
  centrifuge.v_a = 320.;
  centrifuge.height = 1.8;
  centrifuge.diameter = 0.105;
  centrifuge.feed = 13.5 / 1000 / 1000;
  centrifuge.temp = 320.0;
  centrifuge.flow_internal = 2;

  mbmore::CascadeConfig cascade;
  cascade = mbmore::CascadeConfig(
      centrifuge, design_feed_assay, design_product_assay, design_tails_assay,
      design_feed_flow/secpertimestep, max_centrifuges, precision);
  
  
  std::map<int, mbmore::StageConfig>::iterator it;
  for (it = cascade.stgs_config.begin(); it != cascade.stgs_config.end();
       it++) {
    std::cout << "stg " << it->first;
    std::cout << " FA: " << it->second.feed_assay;
    std::cout << " PA: " << it->second.product_assay;
    std::cout << " TA: " << it->second.tail_assay;
    std::cout << " feed_flow: " << it->second.feed_flow;
    std::cout << " cut: " << it->second.cut;
    std::cout << " alpha: " << it->second.alpha;
    std::cout << " beta: " << it->second.beta;
    std::cout << " machine: " << it->second.n_machines;
    std::cout << std::endl;
  }
  std::cout << "Dsign Feed Flow " << cascade.FeedFlow() * secpertimestep  << std::endl;
  return 0;
}

