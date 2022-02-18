//
// Created by tanjunyao7 on 2022/2/18.
//

#ifndef ORB_SLAM3_GARBAGECOLLECTOR_H
#define ORB_SLAM3_GARBAGECOLLECTOR_H
#include "map"
#include "thread"
#include "KeyFrame.h"
#include "MapPoint.h"

namespace ORB_SLAM3
{
    class Atlas;
    class KeyFrameDatabase;

//    class GarbageCollector
//    {
//        GarbageCollector(GarbageCollector& another) = delete;
//        GarbageCollector operator==(GarbageCollector& another) = delete;
//        inline static GarbageCollector& I()
//        {
//            static GarbageCollector instance;
//            return instance;
//        }
//
//        inline void addKeyframe(KeyFrame* kf)
//        {
//            if(!mmKFBuf.count(kf))
//                mmKFBuf[kf] =1;
//            else{
//                mmKFBuf[kf] ++;
//            }
//        }
//
//
//    private:
//        GarbageCollector(){}
//        std::map<KeyFrame*,int> mmKFBuf;
//        std::map<MapPoint*,int> mmMPBuf;
//    };

    template<typename T,int n>
    void PrepareForDeleting(T* data,std::thread::id id)
    {
        static std::map<T*,std::vector<std::thread::id>> dataCount;
        if(!dataCount.count(data))
            dataCount[data] = {id};
        else
        {
            dataCount[data].push_back(id);

            for(auto dc:dataCount)
            {
                std::cout<<dc.first<<std::endl;
                for(auto i: dc.second)
                    std::cout<<i<<" ";
                std::cout<<std::endl;
            }

            if(dataCount[data].size()==n)
            {
//                dataCount.erase(data);
//                delete data;
            }
        }

    }

}

#endif //ORB_SLAM3_GARBAGECOLLECTOR_H
