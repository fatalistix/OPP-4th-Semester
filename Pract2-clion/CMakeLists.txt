cmake_minimum_required(VERSION 3.24)
project(Pract2_clion C)

set(CMAKE_C_STANDARD 17)

add_executable(Pract2_clion parallel_edited_rewrited_buf.c)
set( CMAKE_C_COMPILER mpicxx )

