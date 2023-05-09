#include <gtest/gtest.h>

#include "model.h"

TEST(ModelAlgorithm, RabinKarpAlgorithm) {
  std::filesystem::path path_0{"test/build/HIV-1_AF033819.3.txt"};
  std::filesystem::path path_1{"test/build/test_1.txt"};
  std::filesystem::path path_2{"test/build/output.txt"};

  ASSERT_TRUE(std::filesystem::exists(path_0));
  ASSERT_TRUE(std::filesystem::exists(path_1));

  std::ofstream fi("test/build/output.txt");
  std::cout.rdbuf(fi.rdbuf());

  ModelAlgorithm analyzer;
  analyzer.RabinKarpAlgorithm(path_0, path_1);

  std::cout.rdbuf(std::cout.rdbuf());

  fi.close();

  std::ifstream fo("test/build/output.txt");
  std::string content((std::istreambuf_iterator<char>(fo)),
                      (std::istreambuf_iterator<char>()));
  std::cout << " ";
  ASSERT_EQ(content, "65 9150 9182 \n");
}

TEST(ModelAlgorithm, NWAlgorithm) {
  std::filesystem::path path{"test/build/test_2.txt"};
  ASSERT_TRUE(std::filesystem::exists(path));

  std::ofstream fi("test/build/output.txt");
  std::cout.rdbuf(fi.rdbuf());

  ModelAlgorithm analyzer;
  analyzer.NWAlgorithm(path);

  std::cout.rdbuf(std::cout.rdbuf());

  fi.close();

  std::ifstream fo("test/build/output.txt");
  std::string content((std::istreambuf_iterator<char>(fo)),
                      (std::istreambuf_iterator<char>()));
  std::cout << " ";
  ASSERT_EQ(content,
            "Max score: 10\nGGGCGACACTCCACCATAGA-\n |||||||| |||||||| | "
            "\n-GGCGACAC-CCACCATACAT\n");
}

TEST(ModelAlgorithm, RegexAlgorithm) {
  std::filesystem::path path{"test/build/test_3.txt"};
  ASSERT_TRUE(std::filesystem::exists(path));

  std::ofstream fi("test/build/output.txt");
  std::cout.rdbuf(fi.rdbuf());

  ModelAlgorithm analyzer;
  analyzer.RegexAlgorithm(path);

  std::cout.rdbuf(std::cout.rdbuf());

  fi.close();

  std::ifstream fo("test/build/output.txt");
  std::string content((std::istreambuf_iterator<char>(fo)),
                      (std::istreambuf_iterator<char>()));
  std::cout << " ";
  ASSERT_EQ(content, "True\n");
}

TEST(ModelAlgorithm, KSimilarAlgorithm) {
  std::filesystem::path path{"test/build/test_4.txt"};
  ASSERT_TRUE(std::filesystem::exists(path));

  std::ofstream fi("test/build/output.txt");
  std::cout.rdbuf(fi.rdbuf());

  ModelAlgorithm analyzer;
  analyzer.KSimilarAlgorithm(path);

  std::cout.rdbuf(std::cout.rdbuf());

  fi.close();

  std::ifstream fo("test/build/output.txt");
  std::string content((std::istreambuf_iterator<char>(fo)),
                      (std::istreambuf_iterator<char>()));
  std::cout << " ";
  ASSERT_EQ(content, "3\n");
}

TEST(ModelAlgorithm, WindowAlgorithm) {
  std::filesystem::path path{"test/build/test_5.txt"};
  ASSERT_TRUE(std::filesystem::exists(path));

  std::ofstream fi("test/build/output.txt");
  std::cout.rdbuf(fi.rdbuf());

  ModelAlgorithm analyzer;
  analyzer.WindowAlgorithm(path);

  std::cout.rdbuf(std::cout.rdbuf());

  fi.close();

  std::ifstream fo("test/build/output.txt");
  std::string content((std::istreambuf_iterator<char>(fo)),
                      (std::istreambuf_iterator<char>()));
  std::cout << " ";
  ASSERT_EQ(content, "GACACCCACCATACAT\n");
}