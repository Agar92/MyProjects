project(Test)
cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0012 NEW)

set(name "Test")
project(${name})
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST} "main.cpp")
