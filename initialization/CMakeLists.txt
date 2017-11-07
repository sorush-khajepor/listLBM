file(GLOB Sources "*.f90")
add_library(initialization STATIC ${Sources})
