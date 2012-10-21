set(SRC
    feast.f90
    feast_dense.f90
    f90_functions_wrapper.f90
    )

add_definitions(-Wno-implicit-interface -Wno-unused-dummy-argument
    -Wno-unused-variable -Wno-unused-parameter -cpp)
add_library(feast ${SRC})

install(TARGETS feast
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
    )
