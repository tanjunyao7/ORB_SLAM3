//
// Created by tanjunyao7 on 2022/2/17.
//

#ifndef ORB_SLAM3_MEMORYMONITOR_H
#define ORB_SLAM3_MEMORYMONITOR_H

#include "stdio.h"
#include "sys/types.h"
#include "unistd.h"
#include "string"

void getCurrentProcessMemory(double& vm_usage, double& resident_set)
{
    pid_t proc_id = getpid();
    std::string stat_file_name = "/proc/"+std::to_string(proc_id)+"/stat";

    vm_usage = 0.0;
    resident_set = 0.0;

    ifstream stat_stream(stat_file_name,ios_base::in); //get info from proc directory
    if(!stat_stream.good())
    {
        printf("wrong process id %d \n",proc_id);
        return ;
    }
    //create some variables to get info
    std::string pid, comm, state, ppid, pgrp, session, tty_nr;
    std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    std::string utime, stime, cutime, cstime, priority, nice;
    std::string O, itrealvalue, starttime;
    unsigned long vsize;
    long rss;
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss;
    stat_stream.close();
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured to use 2MB pages

    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

#endif //ORB_SLAM3_MEMORYMONITOR_H
