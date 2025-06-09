// filepath: /home/raksha/gmsh-dealii-interface/include/debug_utils.h
#ifndef DEBUG_UTILS_H
#define DEBUG_UTILS_H

constexpr bool running_in_debug_mode()
{
#ifdef DEBUG
    return true;
#else
    return false;
#endif
}

#endif // DEBUG_UTILS_H