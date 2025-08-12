#include "ModeShape.h"
#include <iostream>
#include <fstream>

int main() {
    ModeShape shape;
    shape.loadFromVTK("no_mem_mode.vtu");
    shape.printVTK() ;
}