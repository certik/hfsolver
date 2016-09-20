set(SRC
    factor.f
    fft235.f
    kernel.f
    zfft3d.f
    )

if(WITH_MPI)
    set(SRC ${SRC}
        pzfft3dv.f
    )
endif()

add_library(ffte ${SRC})

target_compile_options(ffte PRIVATE -std=legacy -Wno-implicit-interface
    -Wno-unused-parameter)

add_library(hfsolver::ffte ALIAS ffte)

install(TARGETS ffte
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
    )
