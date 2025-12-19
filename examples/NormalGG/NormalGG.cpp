/*
 *   This file is part of the OpenPhase (R) software library.
 *  
 *  Copyright (c) 2009-2025 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2025 OpenPhase Solutions GmbH,
 *                Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *     
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.

 *   File created :   2011
 *   Main contributors :  Reza Darvishi Kamachali; Oleg Shchyglo; Hesham Salama
 *
 */

#include "Settings.h"
#include "RunTimeControl.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "Initializations.h"
#include "BoundaryConditions.h"
#include "InterfaceProperties.h"
#include "H5Interface.h"
#include "Tools/TimeInfo.h"
#include "Tools/MicrostructureAnalysis.h"
#include "DrivingForce.h"

using namespace openphase;

int main(int argc, char *argv[])
{
#ifdef MPI_PARALLEL
    int provided = 0;
    OP_MPI_Init_thread(&argc, &argv, OP_MPI_THREAD_FUNNELED, &provided);
    OP_MPI_Comm_rank(OP_MPI_COMM_WORLD, &MPI_RANK);
    OP_MPI_Comm_size(OP_MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif
	std::string InputFile;
    if(argc > 1)
    {
        InputFile = argv[1];
    }
    else
    {
        std::cerr << "No Inputfile provided, trying to use default ProjectInput.opi" << std::endl;
        InputFile = "ProjectInput.opi";
    }
    Settings                    OPSettings;
    OPSettings.ReadInput(InputFile);

    RunTimeControl              RTC(OPSettings, InputFile);
    PhaseField                  Phi(OPSettings, InputFile);
    DoubleObstacle              DO(OPSettings, InputFile);
    InterfaceProperties         IP(OPSettings, InputFile);
    BoundaryConditions          BC(OPSettings, InputFile);
    DrivingForce                dG(OPSettings, InputFile);
    TimeInfo                    Timer(OPSettings, "Execution Time Statistics");
    
    // Initialize HDF5 output
    H5Interface                 H5;
    H5.OpenFile("", "NormalGG_output.h5");
    // WriteSimulationSettings may fail if existing HDF5 structure is incompatible
    // with the current input (HighFive may throw). Guard to avoid aborting run.
    try {
        H5.WriteSimulationSettings(InputFile);
    }
    catch (...) {
        std::cerr << "Warning: could not write simulation settings to HDF5 (continuing)" << std::endl;
    }

    //generating initial grain structure using Voronoi algorithm
    int number_of_grains = 200;
    size_t GrainsPhase = 0;
    Initializations::VoronoiTessellation(Phi, BC, number_of_grains, GrainsPhase);

    

    std::cout << "Entering the Time Loop!!!" << std::endl;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();
        IP.Set(Phi, BC);
        Timer.SetTimeStamp("IP.Set()");
        // prepare/apply any thermodynamic driving force acting on interfaces
        dG.Clear();
        // Compute curvature-driven contribution by default. To use elastic
        // driving force instead, configure @ElasticProperties and call
        // ElasticProperties::CalculateDrivingForce(...) here.
        DO.CalculateCurvatureDrivingForce(Phi, IP, dG);
        dG.Average(Phi, BC);

        // Merge driving force into phase-field increments
        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
        Timer.SetTimeStamp("CalculatePhaseFieldIncrements");
        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("NormalizeIncrements");
        Phi.MergeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("MergeIncrements");

        /// Output to VTK file
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(OPSettings, RTC.tStep);
            // write driving force fields (per-phase averages)
            dG.WriteVTK(OPSettings, Phi, RTC.tStep);
            MicrostructureAnalysis::WriteGrainsStatistics(Phi, RTC.tStep, "NormalGG_output.h5");
            
            // Write HDF5 visualization data (for post-processing)
            std::vector<H5Interface::Field_t> FieldsToWrite;
            FieldsToWrite.push_back(H5Interface::Field_t("PhaseField", 
                [&Phi](int i, int j, int k) -> std::any {
                    return Phi.Fields(i,j,k).get_max().value;
                }));
            FieldsToWrite.push_back(H5Interface::Field_t("GrainIndex", 
                [&Phi](int i, int j, int k) -> std::any {
                    return (double)Phi.Fields(i,j,k).get_max().index;
                }));
            H5.WriteVisualization(RTC.tStep, OPSettings, FieldsToWrite, 1);
            
            // Note: Grain statistics (including HDF5 output) are now written in
            // MicrostructureAnalysis::WriteGrainsStatistics() called above
        }
        /// Output raw data
        if (RTC.WriteRawData())
        {
            Phi.Write(OPSettings, RTC.tStep);
        }
        /// Output to screen
        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message = ConsoleOutput::GetStandard("Interface energy density", I_En);
            ConsoleOutput::WriteTimeStep(RTC, message);
            // print driving force diagnostics and write global averages
            dG.PrintDiagnostics();
            dG.AverageGlobal(Phi, RTC.tStep*RTC.dt);
            Timer.PrintWallClockSummary();
        }
    }
#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize ();
#endif
    return EXIT_SUCCESS;
}
