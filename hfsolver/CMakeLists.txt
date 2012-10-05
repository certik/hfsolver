include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${hfsolver_BINARY_DIR}/src)

cython_add_module_pyx(special)
add_python_library(special special.cpp)
target_link_libraries(special hfsolver)

set(PY_PATH ${PYTHON_INSTALL_PATH}/hfsolver/)
install(TARGETS special
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
