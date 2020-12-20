The script FindNetCDF.cmake was downloaded from
https://github.com/jedbrown/cmake-modules/blob/master/FindNetCDF.cmake
on Dec 20, 2020. I added some hints to the "find_path" and "find_library" commands at the top.

I found another FindNetCDF.cmake file at
https://github.com/Kitware/VTK/blob/master/CMake/FindNetCDF.cmake
That one worked on my laptop and the CI, but not on Draco or Cori.