/*------------------------------------------------------------------------
  Turf: Configurable C++ platform adapter
  Copyright (c) 2016 Jeff Preshing

  Distributed under the Simplified BSD License.
  Original location: https://github.com/preshing/turf

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the LICENSE file for more information.
------------------------------------------------------------------------*/

#include <turf/Core.h>

#if TURF_KERNEL_LINUX

#include <turf/impl/Affinity_Linux.h>
#include <string>
#include <fstream>
#include <stdio.h>

namespace turf {

Affinity_Linux::Affinity_Linux() : m_isAccurate(false), m_numHWThreads(0) {
    std::ifstream f("/proc/cpuinfo");
    if (f.is_open()) {
        CoreInfoCollector collector;
        while (!f.eof()) {
            std::string line;
            std::getline(f, line);
            size_t colon = line.find_first_of(':');
            if (colon != std::string::npos) {
                size_t endLabel = line.find_last_not_of(" \t", colon > 0 ? colon - 1 : 0);
                std::string label = line.substr(0, endLabel != std::string::npos ? endLabel + 1 : 0);
                int value;
                bool isIntegerValue = (sscanf(line.c_str() + colon + 1, " %d", &value) > 0);
                if (isIntegerValue && label == "processor") {
                    collector.flush(*this);
                    collector.logicalProcessor = (s32) value;
                } else if (isIntegerValue && label == "physical id") {
                    collector.coreID.physical = (s32) value;
                } else if (isIntegerValue && label == "core id") {
                    collector.coreID.core = (s32) value;
                }
            }
        }
        collector.flush(*this);
    }
    m_isAccurate = (m_numHWThreads > 0);
    if (!m_isAccurate) {
        m_coreIndexToInfo.resize(1);
        m_coreIndexToInfo[0].hwThreadIndexToLogicalProcessor.push_back(0);
        m_numHWThreads = 1;
    }
}

bool Affinity_Linux::setAffinity(ureg core, ureg hwThread) {
    u32 logicalProcessor = m_coreIndexToInfo[core].hwThreadIndexToLogicalProcessor[hwThread];
    cpu_set_t cpuSet;
    CPU_ZERO(&cpuSet);
    CPU_SET(logicalProcessor, &cpuSet);
    int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpuSet), &cpuSet);
    return (rc == 0);
}

} // namespace turf

#endif // TURF_KERNEL_LINUX
