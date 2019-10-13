#include "testing_utils.hpp"
#include <gtest/gtest.h>

#include <filesystem>

#include "io.hpp"

class TestReceiverFileReading : public testing::Test {
protected:
    std::filesystem::path p = current_source_path(__FILE__) / "data" / "sample_receiverfile";
    std::vector<position_t> expected_result = {
        {5272.000000, 3090.000000, 0.000000}, {5249.959596, 3113.636364, 0.000000},
        {5227.919192, 3137.272727, 0.000000}, {5205.878788, 3160.909091, 0.000000},
        {5183.838384, 3184.545455, 0.000000}, {5161.797980, 3208.181818, 0.000000},
        {5139.757576, 3231.818182, 0.000000}, {5117.717172, 3255.454545, 0.000000},
        {5095.676768, 3279.090909, 0.000000}, {5073.636364, 3302.727273, 0.000000}};
    std::vector<position_t> actual_result = read_receiverfile(p);
};

TEST_F(TestReceiverFileReading, TestIfValuesAreReadCorrectly) {
    // Read values and compare with manually parsed ones
    ASSERT_EQ(actual_result.size(), expected_result.size())
        << "Number of receiver read from receiver file is different than expected.";
    for (auto i = 0; i < expected_result.size(); ++i) {
        EXPECT_EQ(actual_result[i], expected_result[i]) << "Different value for receiver " << i + 1;
    }
}

TEST_F(TestReceiverFileReading, TestIfNonExisitingFileThrows) {
    ASSERT_THROW(read_receiverfile("Idonotexist"), std::runtime_error);
}


class TestSourceFileReading : public testing::Test {
protected:
    std::filesystem::path p = current_source_path(__FILE__) / "data" / "sample_sourcefile";
    std::vector<position_t> expected_result = {{11200, 5600, 10}, {4260, 4407, 10}};
    std::vector<position_t> actual_result = read_sourcefile(p);
};

TEST_F(TestSourceFileReading, TestIfValuesAreReadCorrectly) {
    // Read values and compare with manually parsed ones
    ASSERT_EQ(actual_result.size(), expected_result.size())
        << "Number of sources read from source file different than expected.";
    for (auto i = 0; i < expected_result.size(); ++i) {
        EXPECT_EQ(actual_result[i], expected_result[i]) << "Different value for source " << i;
    }
}

TEST_F(TestSourceFileReading, TestForThrowWhenFileNotExists) {
    ASSERT_THROW(read_sourcefile("Idonotexist"), std::runtime_error);
}