#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string_view>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class DNA_Analyzer {
 public:
  DNA_Analyzer() = default;
  DNA_Analyzer(const DNA_Analyzer& other) = delete;
  DNA_Analyzer(DNA_Analyzer&& other) = delete;
  DNA_Analyzer& operator=(const DNA_Analyzer& other) = delete;
  DNA_Analyzer& operator=(DNA_Analyzer&& other) = delete;
  ~DNA_Analyzer() = default;

  void RabinKarpAlgorithm(const std::filesystem::path& path_1, const std::filesystem::path& path_2);
  void NWAlgorithm(const std::filesystem::path& path);
  void RegexAlgorithm(const std::filesystem::path& path);
  void KSimilarAlgorithm(const std::filesystem::path& path);
  void WindowAlgorithm(const std::filesystem::path& path);

 private:
  const std::size_t rk_mod = 1e9 + 9;
  const std::size_t rk_base = 31;

  std::string ReadFile(const std::filesystem::path& path);
  void CalculateMassHash(const std::string &src, std::vector<std::size_t> &hash, const std::size_t size, const std::vector<std::size_t> &exp);




  struct Score {
    int match, mismatch, gap;
    Score(int match, int mismatch, int gap)
        : match(match), mismatch(mismatch), gap(gap) {}
  };

  int optimalScore(Score score, std::string s1, std::string s2);
  std::pair<int, std::pair<std::string, std::string>> optimalAlignment(Score score, std::string s1, std::string s2);
  bool isMatch(std::string s, std::string p);
  int kSimilarity(std::string s1, std::string s2);
  std::string minWindowSubstring(std::string s, std::string t);
};
