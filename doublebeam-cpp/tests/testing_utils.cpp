#include "testing_utils.hpp"

std::filesystem::path current_source_path(std::string file) {
    std::filesystem::path p(file);
    p = p.parent_path();
    return p;
}