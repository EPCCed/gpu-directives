#include "timer.h"

std::ostream& operator<<(std::ostream & os, const timer & t)
{
    os << t.get_name() << ": ";
    if (t.get_niterations() > 1) 
    {
        os << t.get_time()*1000./t.get_niterations() << " ms/it"  ;
    }
    else
    {
        os << t.get_time() << "s";
    }

    return os;
}