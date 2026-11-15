#pragma once

// Resident memory of the process, for the statistics that the commands of FMSI 2 print: the memory in use
// now, e.g. once the index is loaded or while the queries run, and the peak so far.

#include <cstdint>
#include <cstdio>

#include <sys/resource.h>
#include <unistd.h>

#ifdef __APPLE__
#include <mach/mach.h>
#endif

namespace fmsi2 {

/// Resident memory now, in bytes (0 if the system does not tell).
inline uint64_t resident_bytes() {
#ifdef __APPLE__
    mach_task_basic_info_data_t info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count) != KERN_SUCCESS) return 0;
    return info.resident_size;
#else
    // /proc/self/statm: total and resident pages.
    FILE *statm = fopen("/proc/self/statm", "r");
    if (!statm) return 0;
    unsigned long long total = 0, resident = 0;
    const int read = fscanf(statm, "%llu %llu", &total, &resident);
    fclose(statm);
    return read == 2 ? (uint64_t)resident * (uint64_t)sysconf(_SC_PAGESIZE) : 0;
#endif
}

/// Peak resident memory so far, in bytes.
inline uint64_t peak_resident_bytes() {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#ifdef __APPLE__
    return (uint64_t)usage.ru_maxrss;
#else
    return (uint64_t)usage.ru_maxrss * 1024;
#endif
}

} // namespace fmsi2
