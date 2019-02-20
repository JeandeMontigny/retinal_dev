#ifndef UTIL_NETHODS_
#define UTIL_NETHODS_

namespace bdm {
using namespace std;

  // define my cell creator
  template <typename TCell>
  static void CellCreator(double min, double max, int num_cells, int cellType) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* random = sim->GetRandom();

    for (int i = 0; i < num_cells; i++) {
      double x = random->Uniform(min + 100, max - 100);
      double y = random->Uniform(min + 100, max - 100);
      // RGCL thickness before cell death ~24
      double z = random->Uniform(min + 20, min + 35);
      std::array<double, 3> position = {x, y, z};

      TCell* cell = new TCell(position);
      cell->SetDiameter(random->Uniform(7, 8));  // random diameter
      cell->SetCellType(cellType);
      cell->SetInternalClock(0);
      cell->SetPreviousPosition(position);
      cell->AddBiologyModule(new Substance_secretion_BM());
      cell->AddBiologyModule(new RGC_mosaic_BM());
      cell->AddBiologyModule(new Neurite_creation_BM());
      rm->push_back(cell);
    }
  }  // end CellCreator


  // position exporteur
  inline void position_exporteur(int i) {
    int seed = Simulation::GetActive()->GetRandom()->GetSeed();
    auto* param = Simulation::GetActive()->GetParam();
    ofstream positionFileOn;
    ofstream positionFileOff;
    ofstream positionFileOnOff;
    ofstream positionFileAll;
    stringstream sstmOn;
    stringstream sstmOff;
    stringstream sstmOnOff;
    stringstream sstmAll;
    sstmOn << Concat(param->output_dir_, "/cells_position/on_t") << i << "_seed"
           << seed << ".txt";
    sstmOff << Concat(param->output_dir_, "/cells_position/off_t") << i << "_seed"
            << seed << ".txt";
    sstmOnOff << Concat(param->output_dir_, "/cells_position/onOff_t") << i
              << "_seed" << seed << ".txt";
    sstmAll << Concat(param->output_dir_, "/cells_position/all_t") << i << "_seed"
            << seed << ".txt";

    string fileNameOn = sstmOn.str();
    string fileNameOff = sstmOff.str();
    string fileNameOnOff = sstmOnOff.str();
    string fileNameAll = sstmAll.str();
    positionFileOn.open(fileNameOn);
    positionFileOff.open(fileNameOff);
    positionFileOnOff.open(fileNameOnOff);
    positionFileAll.open(fileNameAll);

    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();

    rm->ApplyOnAllElements([&](SimObject* so) {
      if (auto* soma = so->As<MyCell>()) {
        auto thisCellType = soma->GetCellType();
        auto position = soma->GetPosition();
        if (thisCellType == 0) {
          positionFileOn << position[0] << " " << position[1] << " " << position[2]
                         << "\n";
          positionFileAll << position[0] << " " << position[1] << " " << position[2]
                          << " on\n";
        } else if (thisCellType == 1) {
          positionFileOff << position[0] << " " << position[1] << " " << position[2]
                          << "\n";
          positionFileAll << position[0] << " " << position[1] << " " << position[2]
                          << " off\n";
        } else if (thisCellType == 2) {
          positionFileOnOff << position[0] << " " << position[1] << " "
                            << position[2] << "\n";
          positionFileAll << position[0] << " " << position[1] << " " << position[2]
                          << " onoff\n";
        } else {
          positionFileAll << position[0] << " " << position[1] << " " << position[2]
                          << " nd\n";
        }
      }
    });

    // auto my_cells = rm->template Get<MyCell>();
    // int numberOfCells = my_cells->size();
    //
    // for (int cellNum = 0; cellNum < numberOfCells; cellNum++) {
    //   auto thisCellType = (*my_cells)[cellNum].GetCellType();
    //   auto position = (*my_cells)[cellNum].GetPosition();
    //   if (thisCellType == 0) {
    //     positionFileOn << position[0] << " " << position[1] << " " << position[2]
    //                    << "\n";
    //     positionFileAll << position[0] << " " << position[1] << " " << position[2]
    //                     << " on\n";
    //   } else if (thisCellType == 1) {
    //     positionFileOff << position[0] << " " << position[1] << " " << position[2]
    //                     << "\n";
    //     positionFileAll << position[0] << " " << position[1] << " " << position[2]
    //                     << " off\n";
    //   } else if (thisCellType == 2) {
    //     positionFileOnOff << position[0] << " " << position[1] << " "
    //                       << position[2] << "\n";
    //     positionFileAll << position[0] << " " << position[1] << " " << position[2]
    //                     << " onoff\n";
    //   } else {
    //     positionFileAll << position[0] << " " << position[1] << " " << position[2]
    //                     << " nd\n";
    //   }
    // }

    positionFileOn.close();
    positionFileOff.close();
    positionFileOnOff.close();
    positionFileAll.close();
  }  // end position_exporteur


  template <typename T>
  inline string swc_neurites(T ne, int labelParent,
                             array<double, 3> somaPosition) {
    array<double, 3> nePosition = ne->GetPosition();
    nePosition[0] = nePosition[0] - somaPosition[0];
    nePosition[1] = nePosition[1] - somaPosition[1];
    nePosition[2] = nePosition[2] - somaPosition[2];
    string temps;

    ne->GetMySoma()->IncreaseLabel();
    // set explicitly the value of GetLabel() other wise it is not properly set
    int currentLabel = ne->GetMySoma()->GetLabel();

    // if branching point
    if (ne->GetDaughterRight() != nullptr) {
      // FIXME: segment indice should be 5, no 3. If set here,
      // it's not the actual branching point, but the following segment
      // need to run correction.py to correct file
      temps =
          Concat(temps, "\n", currentLabel, " 3 ", nePosition[0], " ",
                 nePosition[1], " ", nePosition[2], " ", ne->GetDiameter() / 2,
                 " ", labelParent,
                 swc_neurites(ne->GetDaughterRight(), currentLabel, somaPosition))
              .c_str();
      ne->GetMySoma()->IncreaseLabel();
    }
    // if is straigh dendrite
    // need to update currentLabel
    currentLabel = ne->GetMySoma()->GetLabel();
    if (ne->GetDaughterLeft() != nullptr) {
      temps =
          Concat(temps, "\n", currentLabel, " 3 ", nePosition[0], " ",
                 nePosition[1], " ", nePosition[2], " ", ne->GetDiameter() / 2,
                 " ", labelParent,
                 swc_neurites(ne->GetDaughterLeft(), currentLabel, somaPosition))
              .c_str();
    }
    // if ending point
    if (ne->GetDaughterLeft() == nullptr && ne->GetDaughterRight() == nullptr) {
      temps = Concat(temps, "\n", currentLabel, " 6 ", nePosition[0], " ",
                     nePosition[1], " ", nePosition[2], " ",
                     ne->GetDiameter() / 2, " ", labelParent)
                  .c_str();
    }

    return temps;
  }  // end swc_neurites


  inline void morpho_exporteur() {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* param = sim->GetParam();
    int seed = sim->GetRandom()->GetSeed();

    rm->ApplyOnAllElements([&](auto&& so, SoHandle) {
      // FIXME: As<MyCell> ??
      if (auto cell = so->template As<MyNeurite>()) {

        int thisCellType = cell->GetCellType();
        auto cellPosition = cell->GetPosition();
        ofstream swcFile;
        string swcFileName = Concat(param->output_dir_, "/swc_files/cell", cell->GetUid(),
                                    "_type", thisCellType, "_seed", seed, ".swc").c_str();
        swcFile.open(swcFileName);
        cell->SetLabel(1);
        // swcFile << labelSWC_ << " 1 " << cellPosition[0] << " "
        //         << cellPosition[1]  << " " << cellPosition[2] << " "
        //         << cell->GetDiameter()/2 << " -1";
        swcFile << cell->GetLabel() << " 1 0 0 0 " << cell->GetDiameter() / 2
                << " -1";

        for (auto& ne : cell->GetDaughters()) {
          // FIXME: call of swc_neurites() fail
          swcFile << swc_neurites(ne, 1, cellPosition);
        }  // end for neurite in cell
        swcFile.close();
      }
    });  // end for cell in simulation
    std::cout << "swc export done" << std::endl;
  }  // end morpho_exporteur


  // RI computation
  inline double computeRI(vector<array<double, 3>> coordList) {
    vector<double> shortestDistList;
    for (unsigned int i = 0; i < coordList.size();
         i++) {  // for each cell of same type in the simulation
      array<double, 3> cellPosition = coordList[i];

      vector<double> shortestDist;
      // for each other cell of same type in the simulation
      for (unsigned int j = 0; j < coordList.size(); j++) {
        array<double, 3> otherCellPosition = coordList[j];

        // get the distance between those 2 cells (x-y plan only)
        double tempsDistance =
            sqrt(pow(cellPosition[0] - otherCellPosition[0], 2) +
                 pow(cellPosition[1] - otherCellPosition[1], 2));
        // if cell is closer and is not itself
        if (tempsDistance != 0) {
          shortestDist.push_back(tempsDistance);
        }
      }
      // save the shortest distance to neighbour
      shortestDistList.push_back(
        *min_element(shortestDist.begin(), shortestDist.end()));
    }
    // compute mean
    double temps_sum = 0;
    for (unsigned int i = 0; i < shortestDistList.size(); i++) {
      temps_sum += shortestDistList[i];
    }
    double aveShortestDist = temps_sum / (double)shortestDistList.size();
    // compute std
    double temp = 0;
    for (unsigned int i = 0; i < shortestDistList.size(); i++) {
      double val = shortestDistList[i];
      double squrDiffToMean = pow(val - aveShortestDist, 2);
      temp += squrDiffToMean;
    }
    double meanOfDiffs = temp / (double)(shortestDistList.size());
    double std = sqrt(meanOfDiffs);

    // return RI
    return aveShortestDist / std;
  }  // end computeRI


  // get cells coordinate of same cell_type_ to call computeRI
  inline double getRI(int desiredCellType) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    vector<array<double, 3>> coordList;

    rm->ApplyOnAllElements([&](SimObject* so) {
      if (auto* soma = so->As<MyCell>()) {
        auto thisCellType = soma->GetCellType();
        // if cell is of the desired type
        if (thisCellType == desiredCellType) {
          // put cell coord in the list
          coordList.push_back(soma->GetPosition());
        }
      }
    });
    // return RI for desired cell type
    return computeRI(coordList);
  }  // end getRI


  // write file with migration distance of every cells
  inline void exportMigrationDist() {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    auto* param = sim->GetParam();
    vector<array<double, 3>> coordList;
    ofstream migrationDist_outputFile;
    migrationDist_outputFile.open(Concat(param->output_dir_, "/migration_distance.txt"));

    rm->ApplyOnAllElements([&](SimObject* so) {
      if (auto* soma = so->As<MyCell>()) {
        // array<double, 3> positionAtCreation = (*my_cells)[cellNum].GetPreviousPosition();
        // array<double, 3> currentPosition = (*my_cells)[cellNum].GetPosition();
        // double distance = sqrt(pow(currentPosition[0] - positionAtCreation[0], 2) +
        //                        pow(currentPosition[1] - positionAtCreation[1], 2));

        double distance = soma->GetDistanceTravelled();
        migrationDist_outputFile << distance << "\n";
      }
    });

    migrationDist_outputFile.close();
    std::cout << "migration distance export done" << std::endl;
  } // end exportMigrationDist

} // end namespace bdm

#endif
