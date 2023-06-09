#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class ModelAlgorithm {
 public:
  ModelAlgorithm() = default;
  ModelAlgorithm(const ModelAlgorithm& other) = delete;
  ModelAlgorithm(ModelAlgorithm&& other) = delete;
  ModelAlgorithm& operator=(const ModelAlgorithm& other) = delete;
  ModelAlgorithm& operator=(ModelAlgorithm&& other) = delete;
  ~ModelAlgorithm() = default;

  void RabinKarpAlgorithm(const std::filesystem::path& path_1, const std::filesystem::path& path_2);
  void NWAlgorithm(const std::filesystem::path& path);
  void RegexAlgorithm(const std::filesystem::path& path);
  void KSimilarAlgorithm(const std::filesystem::path& path);
  void WindowAlgorithm(const std::filesystem::path& path);

 private:
  const std::array<char, 4> m_ABC = {'A', 'C', 'G', 'T'};
  const std::size_t m_rk_mod = 1e9 + 9;
  const std::size_t m_rk_base = 31;

  std::array<int, 256> m_weights = {'\0'};

  std::string ReadFile(const std::filesystem::path& path);
  void CalculateMassHash(const std::string& src, std::vector<std::size_t>& hash, const std::size_t size, const std::vector<std::size_t>& exp) const;
  std::pair<std::array<int, 3>, std::array<std::string, 2>> GetNWConfig(const std::filesystem::path& path);
  std::vector<std::vector<int>> InitAlighmentMatrix(const std::array<int, 3>& conf, const std::array<std::string, 2>& seq, const std::vector<std::vector<int>>& score_mat);
  std::array<std::string, 3> ComputeAlignedSequences(const std::array<int, 3>& conf, const std::array<std::string, 2>& seq, const std::vector<std::vector<int>>& mat, const std::vector<std::vector<int>>& score_mat) const;
  std::vector<std::vector<int>> InitScoreMatrix(const std::array<int, 3>& conf, const std::size_t rows, const std::size_t cols);
  void PrintNWResult(const std::array<std::string, 3>& alig, const int& score);
  void InitializeWeights();
  char symbol(char c);
};
