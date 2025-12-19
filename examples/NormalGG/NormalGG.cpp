/*
 *   This file is part of the OpenPhase (R) software library.
 *  
 *  Copyright (c) 2009-2025 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2025 OpenPhase Solutions GmbH,
 *                Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *   File created :   2011
 *   Main contributors :  Reza Darvishi Kamachali; Oleg Shchyglo; Hesham Salama
 * 
 * 
 * Modified by :  Yanglang Yuan 
 * 2025-present
 * IBF, RWTH Aachen University, Germany
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
#include <cstdlib>
#include "ConsoleOutput.h"

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
    // Initialize logging as early as possible so that messages from
    // ReadInput and object constructors are captured in the log file.
    ConsoleOutput::InitLogFile("NormalGG.log", VerbosityLevels::Warning);

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
    try {
        H5.WriteSimulationSettings(InputFile);
    }
    catch (...) {
        std::cerr << "Warning: could not write simulation settings to HDF5 (continuing)" << std::endl;
    }

    // determine HDF5 write frequency from Settings (ProjectInput.opi) or default
    int hdf5Freq = OPSettings.HDF5Freq;

    //generating initial grain structure using Voronoi algorithm
    int number_of_grains = 200;
    size_t GrainsPhase = 0;
    Initializations::VoronoiTessellation(Phi, BC, number_of_grains, GrainsPhase);

    

    std::cerr << "Starting..." << std::endl;
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
            
            // Note: Grain statistics and HDF5 visualization are written below
            // under an independent HDF5 frequency (controlled by HDF5_FREQ env).
        }

        // HDF5 output: independent frequency controlled by HDF5_FREQ env var
        if ((RTC.tStep % hdf5Freq) == 0)
        {
            // Write global features to HDF5 (moved to MicrostructureAnalysis for reuse)
            try {
                MicrostructureAnalysis::WriteGlobalFeatures(Phi, DO, H5, RTC, IP);
            } catch (...) {
                std::cerr << "Warning: could not write GlobalFeatures to HDF5 (continuing)" << std::endl;
            }

            // Write driving force fields to HDF5 (per-phase averages)
            try {
                dG.WriteH5(H5, OPSettings, Phi, RTC.tStep);
            } catch (...) {
                std::cerr << "Warning: could not write DrivingForce to HDF5 (continuing)" << std::endl;
            }

            // write grain statistics (HDF5)
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
                        if (OPSettings.WriteDrivingForceH5)
                        {
                            for (size_t a=0; a < Phi.Nphases; ++a)
                                for (size_t b=a; b < Phi.Nphases; ++b)
                                {
                                    std::string name = "dGavg_" + std::to_string(a) + "_" + std::to_string(b);
                                    FieldsToWrite.push_back(H5Interface::Field_t(name,
                                        [&dG,a,b](int i,int j,int k) -> std::any {
                                            return dG.Force(i,j,k).get_average(a,b);
                                        }));
                                }
                        }
            H5.WriteVisualization(RTC.tStep, OPSettings, FieldsToWrite, 1);
        }
        /// Output raw data
        if (RTC.WriteRawData())
        {
            Phi.Write(OPSettings, RTC.tStep);
        }
        /// Output to screen (console) and write diagnostics according to
        /// the configured console output interval. If you want the log to
        /// always contain per-step diagnostics, remove the surrounding
        /// "if (RTC.WriteToScreen())" check.
        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message = ConsoleOutput::GetStandard("Interface energy density", I_En);
            // ConsoleOutput::WriteTimeStep will write to the log and to the
            // console depending on verbosity. We call it only at allowed
            // screen intervals so the log matches the screen interval.
            ConsoleOutput::WriteTimeStep(RTC, message);
            // Driving force diagnostics, global averages and timer summaries use
            // ConsoleOutput internally and therefore will be written to the log
            // only when this block executes.
            dG.PrintDiagnostics();
            dG.AverageGlobal(Phi, RTC.tStep*RTC.dt);
            Timer.PrintWallClockSummary();
        }
    }
#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize ();
#endif
    std::cerr << "Finished" << std::endl;
    ConsoleOutput::CloseLogFile();
    return EXIT_SUCCESS;
}
