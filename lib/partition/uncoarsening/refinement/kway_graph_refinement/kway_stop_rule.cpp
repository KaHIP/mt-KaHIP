#include "kway_stop_rule.h"

#ifdef OUTPUT_GLOBAL_STAT
std::vector<uint32_t> kway_adaptive_stop_rule::m_stat_movements(0);
std::vector<Gain> kway_adaptive_stop_rule::m_stat_gains(0);

std::vector<uint32_t> kway_chernoff_adaptive_stop_rule::m_stat_movements(0);
std::vector<Gain> kway_chernoff_adaptive_stop_rule::m_stat_gains(0);
#endif
