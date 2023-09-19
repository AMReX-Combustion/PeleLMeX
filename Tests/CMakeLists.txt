#=============================================================================
# PeleLMeX Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

if(PELELMEX_ENABLE_FCOMPARE_FOR_TESTS)
  if("${PELELMEX_REFERENCE_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To reference gold files, PELELMEX_REFERENCE_GOLDS_DIRECTORY must be set and exist")
  else()
    if(EXISTS ${PELELMEX_REFERENCE_GOLDS_DIRECTORY})
      set(GOLD_FILES_DIRECTORY ${PELELMEX_REFERENCE_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
      message(STATUS "Test golds directory for fcompare: ${GOLD_FILES_DIRECTORY}")
    else()
      message(FATAL_ERROR "Specified directory for reference gold files does not exist: ${PELELMEX_REFERENCE_GOLDS_DIRECTORY}")
    endif()
  endif()
endif()

if(PELELMEX_SAVE_GOLDS)
  if("${PELELMEX_SAVED_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To save gold files, PELELMEX_SAVED_GOLDS_DIRECTORY must be set and exist")
  else()
    if(EXISTS ${PELELMEX_SAVED_GOLDS_DIRECTORY})
      set(SAVED_GOLDS_DIRECTORY ${PELELMEX_SAVED_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
      message(STATUS "Gold files will be saved to: ${SAVED_GOLDS_DIRECTORY}")
    else()
      message(FATAL_ERROR "Specified directory for saving gold files does not exist: ${PELELMEX_SAVED_GOLDS_DIRECTORY}")
    endif()
  endif()
endif()

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

macro(setup_test)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/Exec/RegTests/${TEST_EXE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR}/tests/${TEST_NAME})
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR}/PeleLMeX-${TEST_EXE_DIR})
    # Gold files should be submodule organized by machine and compiler (these are output during configure)
    set(PLOT_GOLD ${GOLD_FILES_DIRECTORY}/${TEST_EXE_DIR}/tests/${TEST_NAME}/plt00010)
    # Test plot is currently expected to be after 10 steps
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    # Find fcompare
    if(PELELMEX_ENABLE_FCOMPARE_FOR_TESTS)
      set(FCOMPARE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/fcompare)
    endif()
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Gather all files in source directory for test
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*.dat" "${CURRENT_TEST_SOURCE_DIR}/*.py")
    # Copy files to test working directory
    file(COPY ${CURRENT_TEST_SOURCE_DIR}/${TEST_NAME}.inp DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Set some default runtime options for all tests
    set(RUNTIME_OPTIONS "amr.plot_file=plt amr.plot_int=-1 amr.check_int=-1 amrex.the_arena_is_managed=0")
    if(PELELMEX_ENABLE_FPE_TRAP_FOR_TESTS)
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=1 amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1")
    else()
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=0")
    endif()
    if(PELELMEX_ENABLE_MPI)
      if(PELELMEX_ENABLE_CUDA)
        set(PELELMEX_NP 2) # 1 rank per GPU on Eagle
      else()
        set(PELELMEX_NP 4)
      endif()
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PELELMEX_NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(PELELMEX_NP 1)
      unset(MPI_COMMANDS)
    endif()
    # Use fcompare to test diffs in plots against gold files
    if(PELELMEX_ENABLE_FCOMPARE_FOR_TESTS AND (NOT "${TEST_NAME}" MATCHES "hdf5$"))
      if(PELELMEX_ENABLE_CUDA)
        set(FCOMPARE_TOLERANCE "-r 1e-12 --abs_tol 1.0e-12")
      endif()
      set(FCOMPARE_COMMAND "&& ${MPI_COMMANDS} ${FCOMPARE} ${FCOMPARE_TOLERANCE} ${PLOT_TEST} ${PLOT_GOLD}")
    endif()
    if(PELELMEX_SAVE_GOLDS)
      file(MAKE_DIRECTORY ${SAVED_GOLDS_DIRECTORY}/${TEST_EXE_DIR}/${TEST_NAME})
      set(SAVE_GOLDS_COMMAND "&& cp -R ${PLOT_TEST} ${SAVED_GOLDS_DIRECTORY}/${TEST_EXE_DIR}/${TEST_NAME}/")
    endif()
endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME TEST_EXE_DIR)
    setup_test()
    set(RUNTIME_OPTIONS "max_step=10 ${RUNTIME_OPTIONS}")
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp ${RUNTIME_OPTIONS} > ${TEST_NAME}.log ${SAVE_GOLDS_COMMAND} ${FCOMPARE_COMMAND}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 18000 PROCESSORS ${PELELMEX_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "regression" ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_r)

# Regression test with mass conservation verification
function(add_test_rv TEST_NAME TEST_EXE_DIR)
    setup_test()
    set(RUNTIME_OPTIONS "max_step=10 ${RUNTIME_OPTIONS}")
    add_test(${TEST_NAME} sh -c "rm -f datlog && ${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp ${RUNTIME_OPTIONS} > ${TEST_NAME}.log ${SAVE_GOLDS_COMMAND} ${FCOMPARE_COMMAND} && nosetests ${CMAKE_CURRENT_SOURCE_DIR}/test_masscons.py")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 18000 PROCESSORS ${PELELMEX_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "regression;verification" ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_rv)

# Regression tests excluded from CI
function(add_test_re TEST_NAME TEST_EXE_DIR)
    add_test_r(${TEST_NAME} ${TEST_EXE_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;no-ci")
endfunction(add_test_re)

# Regression tests expected to fail
function(add_test_rf TEST_NAME TEST_EXE_DIR)
    add_test_r(${TEST_NAME} ${TEST_EXE_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES WILL_FAIL TRUE)
endfunction(add_test_rf)

# Verification test with 1 resolution
function(add_test_v1 TEST_NAME TEST_EXE_DIR)
    setup_test()
    set(RUN_COMMAND "rm -f mmslog datlog && ${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp")
    add_test(${TEST_NAME} sh -c "${RUN_COMMAND} ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && nosetests ${TEST_NAME}.py")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 18000 PROCESSORS ${PELELMEX_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "verification" ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_v1)

# Verification test with multiple resolutions (each test runs on maximum number of processes on node)
function(add_test_v2 TEST_NAME TEST_EXE_DIR LIST_OF_GRID_SIZES)
    setup_test()
    unset(MASTER_RUN_COMMAND)
    # Get last item in resolution list so we can find out when we are on the last item in our loop
    list(GET LIST_OF_GRID_SIZES -1 LAST_GRID_SIZE_IN_LIST)
    # Create the commands to run for each resolution
    foreach(GRID_SIZE IN LISTS LIST_OF_GRID_SIZES)
      file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE})
      file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*.dat" "${CURRENT_TEST_SOURCE_DIR}/*.py")
      file(COPY ${CURRENT_TEST_SOURCE_DIR}/${TEST_NAME}.inp DESTINATION "${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/")
      file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/")
      if(${PELELMEX_DIM} EQUAL 3)
        set(NCELLS "${GRID_SIZE} ${GRID_SIZE} ${GRID_SIZE}")
      elseif(${PELELMEX_DIM} EQUAL 2)
        set(NCELLS "${GRID_SIZE} ${GRID_SIZE}")
      elseif(${PELELMEX_DIM} EQUAL 1)
        set(NCELLS "${GRID_SIZE}")
      endif()
      set(DELETE_PREVIOUS_FILES_COMMAND "rm -f mmslog datlog")
      set(RUN_COMMAND_${GRID_SIZE} "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/${TEST_NAME}.inp")
      set(RUNTIME_OPTIONS_${GRID_SIZE} "${RUNTIME_OPTIONS} amr.n_cell=${NCELLS}")
      string(APPEND MASTER_RUN_COMMAND "cd ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}")
      string(APPEND MASTER_RUN_COMMAND " && ")
      string(APPEND MASTER_RUN_COMMAND "${DELETE_PREVIOUS_FILES_COMMAND}")
      string(APPEND MASTER_RUN_COMMAND " && ")
      string(APPEND MASTER_RUN_COMMAND "${RUN_COMMAND_${GRID_SIZE}} ${RUNTIME_OPTIONS_${GRID_SIZE}} > ${TEST_NAME}-${GRID_SIZE}.log")
      # Add another " && " unless we are on the last resolution in the list
      if(NOT ${GRID_SIZE} EQUAL ${LAST_GRID_SIZE_IN_LIST})
        string(APPEND MASTER_RUN_COMMAND " && ")
      endif()
    endforeach()
    set(IMAGES_TO_UPLOAD ${CURRENT_TEST_BINARY_DIR}/p_error.png ${CURRENT_TEST_BINARY_DIR}/rho_error.png ${CURRENT_TEST_BINARY_DIR}/u_error.png)
    if(${PELELMEX_DIM} EQUAL 3)
      list(APPEND IMAGES_TO_UPLOAD ${CURRENT_TEST_BINARY_DIR}/v_error.png ${CURRENT_TEST_BINARY_DIR}/w_error.png)
    elseif(${PELELMEX_DIM} EQUAL 2)
      list(APPEND IMAGES_TO_UPLOAD ${CURRENT_TEST_BINARY_DIR}/v_error.png)
    endif()
    add_test(${TEST_NAME} sh -c "${MASTER_RUN_COMMAND} && cd ${CURRENT_TEST_BINARY_DIR} && nosetests ${CMAKE_CURRENT_SOURCE_DIR}/test_second_order.py")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 18000 PROCESSORS ${PELELMEX_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}" LABELS "verification;no-ci" ATTACHED_FILES "${IMAGES_TO_UPLOAD}")
endfunction(add_test_v2)

# Standard unit test
function(add_test_u TEST_NAME)
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/Exec/UnitTests/PeleLMeX-UnitTests)
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/Exec/UnitTests/${TEST_EXE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/Exec/UnitTests/${TEST_EXE_DIR})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    set(PELELMEX_NP 1)
    if(PELELMEX_ENABLE_MPI)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PELELMEX_NP} ${MPIEXEC_PREFLAGS}")
    else()
      unset(MPI_COMMANDS)
    endif()
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 1800 PROCESSORS ${PELELMEX_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "unit")
endfunction(add_test_u)

function(add_test_spray TEST_EXE_DIR)
    set(TEST_NAME ${TEST_EXE_DIR})
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/Exec/RegTests/${TEST_EXE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR})
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR}/PeleLMeX-${TEST_EXE_DIR})
    if(PELELMEX_ENABLE_MPI)
      if(PELELMEX_ENABLE_CUDA)
        set(PELELMEX_NP 2) # 1 rank per GPU on Eagle
      else()
        set(PELELMEX_NP 4)
      endif()
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PELELMEX_NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(PELELMEX_NP 1)
      unset(MPI_COMMANDS)
    endif()
    set(FCOMPARE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/fcompare)
    set(MULTIRUN_FLAGS "--test_dir ${CURRENT_TEST_BINARY_DIR} --run_cmd '${MPI_COMMANDS}' --pele_exec PeleLMeX-${TEST_EXE_DIR} --fcomp_exec ${FCOMPARE} --input_file ${TEST_NAME}.inp")
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*.dat" "${CURRENT_TEST_SOURCE_DIR}/*.py")
    # Copy files to test working directory
    file(COPY ${CURRENT_TEST_SOURCE_DIR}/${TEST_NAME}.inp DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    add_test(${TEST_NAME} sh -c "python3 ${CURRENT_TEST_BINARY_DIR}/multiRuns.py ${MULTIRUN_FLAGS}")
endfunction(add_test_spray)

#=============================================================================
# Regression tests
#=============================================================================

# Run in CI
add_test_r(flamesheet-drm19-${PELELMEX_DIM}d FlameSheet)
add_test_re(covo-${PELELMEX_DIM}d PeriodicCases)
add_test_re(hotbubble-${PELELMEX_DIM}d HotBubble)
add_test_re(enclosedflame-${PELELMEX_DIM}d EnclosedFlame)
if(PELELMEX_DIM EQUAL 3)
  add_test_re(hit-${PELELMEX_DIM}d HITDecay)
  add_test_re(hit-les-${PELELMEX_DIM}d HITDecay)
endif()
