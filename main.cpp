#include <iostream>
#include "./googletest/googletest/include/gtest/gtest.h"

#include <fstream>
#include <filesystem>

int main( int argc, char** argv ) {
    testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
    return 0;

    std::cout << "Hello, World!" << std::endl;
    return 0;
}
