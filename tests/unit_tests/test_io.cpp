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
    ASSERT_EQ(read_timesteps(p), t);
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
    std::vector<Receiver> expected_result = {
        {5272.000000_meter, 3090.000000_meter, 0.000000_meter, 1},
        {5249.959596_meter, 3113.636364_meter, 0.000000_meter, 2},
        {5227.919192_meter, 3137.272727_meter, 0.000000_meter, 3},
        {5205.878788_meter, 3160.909091_meter, 0.000000_meter, 4},
        {5183.838384_meter, 3184.545455_meter, 0.000000_meter, 5},
        {5161.797980_meter, 3208.181818_meter, 0.000000_meter, 6},
        {5139.757576_meter, 3231.818182_meter, 0.000000_meter, 7},
        {5117.717172_meter, 3255.454545_meter, 0.000000_meter, 8},
        {5095.676768_meter, 3279.090909_meter, 0.000000_meter, 9},
        {5073.636364_meter, 3302.727273_meter, 0.000000_meter, 10}};
    std::vector<Receiver> actual_result = read_receiverfile(p);
};

TEST_F(TestReceiverFileReading, TestIfValuesAreReadCorrectly) {
    // Read values and compare with manually parsed ones
    ASSERT_EQ(actual_result.size(), expected_result.size())
        << "Number of receiver read from receiver file is different than expected.";
    std::cout << std::setprecision(17);
    for (auto i = 0; i < expected_result.size(); ++i) {
        EXPECT_EQ(actual_result[i], expected_result[i])
            << "Different value for receiver " << i + 1 << ": " << actual_result[i] << " vs "
            << expected_result[i];
    }
}

TEST_F(TestReceiverFileReading, TestIfNonExisitingFileThrows) {
    ASSERT_THROW(read_receiverfile("Idonotexist"), std::runtime_error);
}


class TestSourceFileReading : public testing::Test {
protected:
    std::filesystem::path p = current_source_path(__FILE__) / "data" / "sample_sourcefile";
    std::vector<Source> expected_result = {{11200_meter, 5600_meter, 10_meter, 1},
                                           {4260_meter, 4407_meter, 10_meter, 2}};
    std::vector<Source> actual_result = read_sourcefile(p);
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

TEST_F(TestSourceFileReading, TestIfMissingNumberOfSourcesThrows) {
    p.replace_filename("sample_sourcefile_missing_nsrc");
    ASSERT_THROW(read_sourcefile(p), std::runtime_error);
}

TEST_F(TestSourceFileReading, TestIfMissingSourceNumbersAreRecognized) {
    p.replace_filename("sample_sourcefile_missing_sourcenumbers");
    ASSERT_THROW(read_sourcefile(p), std::runtime_error);
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