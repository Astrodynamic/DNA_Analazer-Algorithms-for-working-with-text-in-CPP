#pragma once

#include <fstream>
#include <iostream>
#include <filesystem>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class DNA_Analyzer {
 public:
  DNA_Analyzer() = default;
  DNA_Analyzer(const DNA_Analyzer &other) = delete;
  DNA_Analyzer(DNA_Analyzer &&other) = delete;
  DNA_Analyzer &operator=(const DNA_Analyzer &other) = delete;
  DNA_Analyzer &operator=(DNA_Analyzer &&other) = delete;
  ~DNA_Analyzer() = default;

  void RabinKarpAlgorithm(const std::filesystem::path& path_1, const std::filesystem::path& path_2);
  void NWAlgorithm(const std::filesystem::path& path);
  void RegexAlgorithm(const std::filesystem::path& path);
  void KSimilarAlgorithm(const std::filesystem::path& path);
  void WindowAlgorithm(const std::filesystem::path& path);

 private:

 struct Score {
  int match, mismatch, gap;
  Score(int _match, int _mismatch, int _gap) : match(_match), mismatch(_mismatch), gap(_gap) {}
};

  int optimalScore(Score score, std::string s1, std::string s2);
  std::pair<int, std::pair<std::string, std::string>> optimalAlignment(Score score, std::string s1, std::string s2);
  bool isMatch(std::string s, std::string p);
  int kSimilarity(std::string s1, std::string s2);
  std::string minWindowSubstring(std::string s, std::string t);
};
