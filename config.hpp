#ifndef VTKMESH_COMMON_H__
#define VTKMESH_COMMON_H__

#ifndef VTKMESH_COMPILED_LIB
    #define VTKMESH_HEADER_ONLY
#endif

#ifdef VTKMESH_HEADER_ONLY
    #define VTKMESH_INLINE inline
#else
    #define VTKMESH_INLINE
#endif

#if defined(_WIN32) && defined(VTKMESH_COMPILED_LIB)
    #define VTKMESH_API __declspec(dllexport)
#else
    #define VTKMESH_API
#endif

#endif
