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
#include <optional>
#include <array>

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
  const std::array<char, 4> m_ABC = {'A', 'C', 'G', 'T'};
  const std::size_t m_rk_mod = 1e9 + 9;
  const std::size_t m_rk_base = 31;

  std::array<int, 256> m_weights = {};

  std::string ReadFile(const std::filesystem::path& path);
  void CalculateMassHash(const std::string &src, std::vector<std::size_t> &hash, const std::size_t size, const std::vector<std::size_t> &exp);
  std::pair<std::array<int, 3>, std::array<std::string, 2>> GetNWConfig(const std::filesystem::path& path);
  std::vector<std::vector<int> > InitAlighmentMatrix(const std::array<int, 3> &conf, const std::array<std::string, 2> &seq, const std::vector<std::vector<int> > &score_mat);
  std::array<std::string, 3> ComputeAlignedSequences(const std::array<int, 3> &conf, const std::array<std::string, 2> &seq, const std::vector<std::vector<int> > &mat, const std::vector<std::vector<int> > &score_mat);
  std::vector<std::vector<int> > InitScoreMatrix(const std::array<int, 3> &conf, const std::size_t rows, const std::size_t cols);
  void PrintNWResult(const std::array<std::string, 3> &alig, const int &score);
  void InitializeWeights();


  char symbol(char c);
  int kSimilarity(std::string s1, std::string s2);
  std::string minWindowSubstring(std::string s, std::string t);
};
