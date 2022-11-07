#include <iostream>
#include "hello.h"

std::vector<int> hello(){
    std::cout << "HELLO" << std::endl;
    std::vector<int> a(3, 0);
    return a;
}