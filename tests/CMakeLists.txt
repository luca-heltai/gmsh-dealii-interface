add_definitions(-DDEBUG)  
INCLUDE_DIRECTORIES(${D2K_INCLUDE_DIRS})
# Automatically pick up other tests
DEAL_II_PICKUP_TESTS()