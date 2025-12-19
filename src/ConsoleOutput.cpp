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

 *   File created :   2013
 *   Main contributors :   Philipp Engels; Raphael Schiedung
 *
 */

#include "ConsoleOutput.h"
#include "RunTimeControl.h"
#include "Includes.h"
#include "BuildInfo.h"

namespace openphase
{

using namespace std;

void ConsoleOutput::InitLogFile(const std::string& filename, VerbosityLevels consoleVerbosity)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK != 0) return;
#endif
    if (filename.empty())
    {
        LogEnabled = false;
        return;
    }
    LogFile.open(filename, std::ios::out | std::ios::trunc);
    if (LogFile.is_open())
    {
        LogEnabled = true;
        OutputVerbosity = consoleVerbosity;
        // redirect std::cout to the log file so any direct cout goes into log
        CoutBuf = std::cout.rdbuf();
        std::cout.rdbuf(LogFile.rdbuf());
        LogFile << "--- OpenPhase log started: " << get_time() << " ---" << std::endl;
    }
    else
    {
        LogEnabled = false;
    }
}

void ConsoleOutput::CloseLogFile()
{
#ifdef MPI_PARALLEL
    if (MPI_RANK != 0) return;
#endif
    if (LogEnabled)
    {
        LogFile << "--- OpenPhase log finished: " << get_time() << " ---" << std::endl;
        // restore cout rdbuf
        if (CoutBuf)
        {
            std::cout.rdbuf(CoutBuf);
            CoutBuf = nullptr;
        }
        LogFile.close();
        LogEnabled = false;
    }
}

void ConsoleOutput::WriteToLog(const std::string& message)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK != 0) return;
#endif
    if (LogEnabled)
    {
        LogFile << message;
        // Ensure flushed so logs are available during long runs
        LogFile.flush();
    }
}

std::string ConsoleOutput::get_time()
{
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string s(30, '\0');
    size_t size =  std::strftime(&s[0], s.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    s.resize(size);
    return s;
}

void ConsoleOutput::WriteTimeStep(const RunTimeControl& RTC, const string Message,
        size_t ColumnWidth)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(true)
    {
        std::ostringstream oss;
        oss << setfill(' ') << setw(ColumnWidth)  << left << "Time step" << ": " << to_string(RTC.TimeStep) + "/" + to_string(RTC.MaxTimeStep) << "\n";
        oss << setfill(' ') << setw(ColumnWidth)  << left << "Simulation time" << ": " << RTC.SimulationTime << "\n";
        oss << setfill(' ') << setw(ColumnWidth)  << left << "Wall clock time" << ": " << get_time() << "\n";
        if (Message != "")
        {
            oss << Message;
        }
        // Write to log regardless of console verbosity
        WriteToLog(string("===============================\n") + oss.str() + string("===============================\n"));
        // Print to console only if verbosity level permits
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            WriteLine("=");
            cout << oss.str();
            if (Message != "") { WriteLine("-"); cout << Message; }
            WriteLine("=");
        }
    }
}

void ConsoleOutput::WriteTimeStep(const int tStep, const int nSteps, const string Message,
        size_t ColumnWidth)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(true)
    {
        std::ostringstream oss;
        oss << setfill(' ') << setw(ColumnWidth)  << left << "TimeStep" << " " << to_string(tStep) + "/" + to_string(nSteps) << "\n";
        oss << setfill(' ') << setw(ColumnWidth)  << left << "Time" << " " << get_time() << "\n";
        if (Message != "")
        {
            oss << Message << "\n";
        }
        WriteToLog(string("-------------------------------\n") + oss.str() + string("-------------------------------\n"));
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            WriteLine("_");
            cout << oss.str();
            WriteLine("_");
        }
    }
}

void ConsoleOutput::WriteTimeStep(const int tScreenWrite, const int tStep,
        const int nSteps, const string Message, size_t ColumnWidth)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    // Always write time step to log (if enabled). Console output depends on verbosity
    if (!(tStep%tScreenWrite))
    {
        WriteTimeStep(tStep, nSteps, Message, ColumnWidth);
    }
}

void ConsoleOutput::WriteLine(const string LineType)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        std::string line(LineLength, LineType[0]);
        line += "\n";
        WriteToLog(line);
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            cout << line;
        }
    }
}

void ConsoleOutput::WriteLineInsert(const string Insert, const string LineType)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        std::ostringstream oss;
        oss << LineType[0] << LineType[0] << "< " << Insert << " >" << setfill(LineType[0]) << setw(std::max<int>(0,LineLength-Insert.size()-6)) << "" << "\n";
        WriteToLog(oss.str());
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            cout << oss.str();
        }
    }
}

void ConsoleOutput::WriteBlankLine(void)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        WriteToLog("\n");
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            cout << "\n";
        }
    }
}

void ConsoleOutput::WriteSimple(const string Message)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        std::string out = Message + "\n";
        WriteToLog(out);
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            cout << out;
        }
    }
}

void ConsoleOutput::WriteCoordinate(const int x, const int y, const int z, const double dx)
{
    WriteStandard("Point",iVector3{x,y,x});
    WriteStandard("Coordinate",dVector3{x*dx,y*dx,x*dx});
}

void ConsoleOutput::WriteWithinMethod(const string Message, const string Instance,
                                      const string Method)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        std::ostringstream oss;
        string thisInstance = Instance;
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        oss << setfill(' ') << setw(10)  << left << thisInstance << "\n" << Message << "\n";
        oss << string(LineLength, '-') << "\n";
        WriteToLog(oss.str());
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            cout << oss.str();
        }
    }
}

void ConsoleOutput::WriteWarning(const string Message, const string Instance,
                        const string Method)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity >= VerbosityLevels::Warning)
    {
        string thisInstance = Instance;
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        std::ostringstream oss;
        oss << setfill('~') << setw(LineLength) << "" << "\n";
        oss << setfill(' ') << setw(10)  << left << "Warning: "  << thisInstance << "\n"
            << "          " << Message      << "\n";
        oss << setfill('~') << setw(LineLength) << "" << "\n";
        // write to log first
        WriteToLog(oss.str());
        // then to console (cerr)
        cerr << oss.str();
    }
}

void ConsoleOutput::WriteExit(const string Message, const string Instance, const string Method)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        string thisInstance = Instance;
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        std::ostringstream oss;
        oss << "\n";
        oss << setfill('*') << setw(LineLength) << "" << "\n";
        oss << setfill(' ') << setw(10) << left << "Calculation terminated!" << "\n"
            << setw(10) << left << "Instance:" << thisInstance << "\n"
            << setw(10) << left << "Reason: " << Message << "\n"
            << setw(10) << left << "Time: " << get_time() << "\n";
        oss << setfill('*') << setw(LineLength) << "" << "\n";
        oss << "\n";
        WriteToLog(oss.str());
        cerr << oss.str();
    }
}

void ConsoleOutput::WriteStartScreen(void)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(true)
    {
        std::ostringstream oss;
        oss << "\n";
        oss << string(LineLength, '>') << "\n";
        oss << "  OpenPhase\n\n"
            << "  Copyright (c) Ruhr-Universitaet Bochum, Universitaetsstrasse 150, 44801 Bochum, Germany\n"
            << "            and OpenPhase Solutions GmbH, Universitaetsstrasse 136, 44799 Bochum, Germany.\n"
            << "  All rights reserved.\n";
        oss << string(LineLength, '<') << "\n";
        oss << "\n";
        oss << "Build time: " << BUILD_TIME << "\n";
        oss << "Git commit SHA: " << GIT_COMMIT_SHA << "\n";
        oss << std::endl;
        WriteToLog(oss.str());
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            cout << oss.str();
        }
    }
}

void ConsoleOutput::PressEnterToContinue(void)
{
    cout << "Press ENTER to continue... " << flush;
    cin.ignore( numeric_limits <streamsize> ::max(), '\n' );
}

void ConsoleOutput::StartProgressIndicator(std::string IndicatorTitle)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(true)
    {
        std::ostringstream oss;
        oss << IndicatorTitle << ":\n";
        oss << "0------------------25------------------50------------------75---------------100%" << std::endl;
        WriteToLog(oss.str());
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            cout << oss.str();
        }
    }
}

void ConsoleOutput::AdvanceProgressIndicator(double start_value, double end_value, double current_value, int& pos)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(true)
    {
        int current_value_normalized = 80.0*fabs((current_value - start_value)/(end_value - start_value));
        std::string out;
        for(int n = pos; n < current_value_normalized; n++)
        {
            pos++;
            out.push_back('#');
        }
        if(!out.empty())
        {
            WriteToLog(out);
        }
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            cout << out << flush;
        }
    }
}

void ConsoleOutput::EndProgressIndicator(void)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(true)
    {
        WriteToLog("\nDone!\n");
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            cout << endl;
            ConsoleOutput::WriteSimple("Done!");
        }
    }
}

}// namespace openphase
