#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "CascadeConfig.h"
#include "CentrifugeConfig.h"
#include "StageConfig.h"

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cout << "Bad Number of Args, need 3 and only 3!" << std::endl;
    return 1;
  }

  int secpertimestep = 2629846;

  double design_feed_assay = std::atof(argv[1]) * 0.01;
  double design_product_assay = std::atof(argv[2]) * 0.01;
  double design_tails_assay = std::atof(argv[3]) * 0.01;

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
      design_feed_flow / secpertimestep, max_centrifuges, precision);







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
    std::cout << " flow per machine: " << it->second.n_machines* it->second.centrifuge.feed;
    std::cout << " ratio: " << it->second.n_machines* it->second.centrifuge.feed / it->second.feed_flow;
    std::cout << std::endl;
  }
  std::cout << "Dsign Feed Flow " << cascade.FeedFlow() * secpertimestep
            << std::endl;








  double assay = 0.0035;
  double stop_assay = 0.95;
  double assay_increment = 0.001;

  std::ofstream model_0;
  model_0.open("model_0.dat");
  std::ofstream model_1;
  model_1.open("model_1.dat");
  std::ofstream model_2;
  model_2.open("model_2.dat");

  while (assay < stop_assay) {
    mbmore::CascadeConfig cascade_tmp =
        cascade.ModelMissUsedCascade(assay, 0, precision);
    double product_assay =
        cascade_tmp.stgs_config.rbegin()->second.product_assay;
    double tail_assay = cascade_tmp.stgs_config.begin()->second.tail_assay;
    double feed_flow = cascade_tmp.FeedFlow();
    double product_flow = cascade_tmp.ProductFlow();
    double tail_flow = cascade_tmp.TailFlow();
    model_0 << assay << " " << product_assay << " " << tail_assay << " "
            << feed_flow << " " << product_flow << " " << tail_flow
            << std::endl;

    cascade_tmp = cascade.ModelMissUsedCascade(assay, 1, precision);
    product_assay = cascade_tmp.stgs_config.rbegin()->second.product_assay;
    tail_assay = cascade_tmp.stgs_config.begin()->second.tail_assay;
    feed_flow = cascade_tmp.FeedFlow();
    product_flow = cascade_tmp.ProductFlow();
    tail_flow = cascade_tmp.TailFlow();
    model_1 << assay << " " << product_assay << " " << tail_assay << " "
            << feed_flow << " " << product_flow << " " << tail_flow
            << std::endl;

    cascade_tmp = cascade.ModelMissUsedCascade(assay, 2, precision);
    product_assay = cascade_tmp.stgs_config.rbegin()->second.product_assay;
    tail_assay = cascade_tmp.stgs_config.begin()->second.tail_assay;
    feed_flow = cascade_tmp.FeedFlow();
    product_flow = cascade_tmp.ProductFlow();
    tail_flow = cascade_tmp.TailFlow();
    model_2 << assay << " " << product_assay << " " << tail_assay << " "
            << feed_flow << " " << product_flow << " " << tail_flow
            << std::endl;

    assay += assay_increment;
  }
  return 0;
}
