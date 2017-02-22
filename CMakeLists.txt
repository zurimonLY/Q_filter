WRAP_GGO(ds3zipmlsart_GGO_C ds3zipmlsart.ggo ../rtkinputprojections_section.ggo ../rtk3Doutputimage_section.ggo ${RTK_BINARY_DIR}/rtkVersion.ggo)
add_executable(ds3zipmlsart ds3zipmlsart.cxx ${ds3zipmlsart_GGO_C})
target_link_libraries(ds3zipmlsart RTK)
add_definitions(-std=c++11)
if (RTK_USE_CUDA)
  target_link_libraries(ds3zipmlsart rtkcuda)
endif ()

# Installation code
if(NOT RTK_INSTALL_NO_EXECUTABLES)
  foreach(EXE_NAME ds3zipmlsart) 
    install(TARGETS ${EXE_NAME}
      RUNTIME DESTINATION ${RTK_INSTALL_RUNTIME_DIR} COMPONENT Runtime
      LIBRARY DESTINATION ${RTK_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
      ARCHIVE DESTINATION ${RTK_INSTALL_ARCHIVE_DIR} COMPONENT Development)
  endforeach() 
endif()

