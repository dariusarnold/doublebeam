#include "testing_utils.hpp"
#include <gtest/gtest.h>

#include <filesystem>

#include "io.hpp"

class TestSeismogramLoading : public testing::Test {
protected:
    std::filesystem::path p = current_source_path(__FILE__) / "data" / "receiver_001.txt";
    std::vector<double> t{0.000, 0.004, 0.008, 0.012, 0.016, 0.020,
                          0.024, 0.028, 0.032, 0.036, 0.040};
    std::vector<double> x{-6.086347e-27, 4.000257e-28,  6.390645e-27,  5.060195e-27,
                          4.030744e-28,  -1.626664e-27, -7.453548e-28, 9.324629e-28,
                          1.507141e-27,  -1.028619e-27, -5.575252e-27};
};

TEST_F(TestSeismogramLoading, TestIfTimeIsLoadedCorrectly) {
    ASSERT_EQ( read_timesteps(p), t);
}

TEST_F(TestSeismogramLoading, TestIfAmplitudeIsReadCorrectly) {
    ASSERT_EQ(read_amplitude(p), x);
}

TEST_F(TestSeismogramLoading, TestNonExistingFileThrows_Amplitude) {
    ASSERT_THROW(read_timesteps("Idonotexist"), std::runtime_error);
}

TEST_F(TestSeismogramLoading, TestNonExistingFileThrows_Timesteps) {
    ASSERT_THROW(read_timesteps("Idonotexist"), std::runtime_error);
}


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


class TestVectorBinaryIO : public testing::Test {
protected:
    std::vector<double> data{0.5, 1 / 3, 42, 0, -10.1};
    std::filesystem::path p = current_source_path(__FILE__) / "data" / "vector.bin";
};

TEST_F(TestVectorBinaryIO, TestStoringToDisk) {
    if (std::filesystem::exists(p)) {
        // clean up old data
        std::filesystem::remove(p);
    }
    ASSERT_FALSE(std::filesystem::exists(p))
        << "Failed to clean up for unit test (Failed to delete file " << p << ").";
    ASSERT_NO_THROW(save_binary(data, p)) << "Failed to save vector to file " << p << ".";
    ASSERT_TRUE(std::filesystem::exists(p));
}

TEST_F(TestVectorBinaryIO, TestReadingFromDisk) {
    if (std::filesystem::exists(p)) {
        std::filesystem::remove(p);
    }
    ASSERT_FALSE(std::filesystem::exists(p))
        << "Failed to clean up for unit test (Failed to delete file " << p << ").";
    ASSERT_NO_THROW(save_binary(data, p))
        << "Failed to set up test environment (Couldn't save vector file " << p << ").";
    std::vector<double> v;
    ASSERT_NO_THROW(v = load_binary<double>(p)) << "Couldn't load binary vector from file.";
    ASSERT_EQ(v.size(), data.size());
    ASSERT_EQ(v, data);
}