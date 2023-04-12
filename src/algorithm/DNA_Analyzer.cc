#include "DNA_Analyzer.h"

void DNA_Analyzer::RabinKarpAlgorithm(const std::filesystem::path& path_1, const std::filesystem::path& path_2) {
  std::string text = std::move(ReadFile(path_1));
  std::string samp = std::move(ReadFile(path_2));

  const std::size_t text_size = text.size();
  const std::size_t samp_size = samp.size();

  std::vector<std::size_t> exp(std::max(text_size, samp_size), 1);
  std::transform(exp.begin(), exp.end() - 1, exp.begin() + 1, [this](std::size_t& item) { return (item * m_rk_base) % m_rk_mod; });

  std::vector<std::size_t> hash_text(text_size + 1, 0);
  CalculateMassHash(text, hash_text, text_size, exp);

  std::vector<std::size_t> hash_samp(samp_size + 1, 0);
  CalculateMassHash(samp, hash_samp, samp_size, exp);

  for (std::size_t idx = 0; idx + samp_size - 1 < text_size; ++idx) {
    std::size_t hash_curr = (hash_text[idx + samp_size] + m_rk_mod - hash_text[idx]) % m_rk_mod;
    if ((hash_curr == (hash_samp[samp_size] * exp[idx] % m_rk_mod)) && text.substr(idx, samp_size) == samp) {
      std::cout << idx << " ";
    }
  }
  std::cout << std::endl;
}

void DNA_Analyzer::CalculateMassHash(const std::string& src, std::vector<std::size_t>& hash, const std::size_t size, const std::vector<std::size_t>& exp) {
  for (std::size_t i = 0; i < size; ++i) {
    hash[i + 1] = (hash[i] + src[i] * exp[i]) % m_rk_mod;
  }
}

std::string DNA_Analyzer::ReadFile(const std::filesystem::path& path) {
  std::ifstream file(path);
  return std::string((std::istreambuf_iterator<char>(file)), {});
}

void DNA_Analyzer::NWAlgorithm(const std::filesystem::path& path) {
  auto [conf, seq] = GetNWConfig(path);

  InitializeWeights();
  auto score = InitScoreMatrix(conf, 4, 4);
  auto mat = InitAlighmentMatrix(conf, seq, score);
  auto alig = ComputeAlignedSequences(conf, seq, mat, score);
  PrintNWResult(alig, mat[seq[0].length()][seq[1].length()]);
}

std::pair<std::array<int, 3>, std::array<std::string, 2>> DNA_Analyzer::GetNWConfig(const std::filesystem::path& path) {
  std::ifstream file(path);

  std::array<int, 3> conf;
  std::array<std::string, 2> seq;

  file >> conf[0] >> conf[1] >> conf[2];
  file >> seq[0] >> seq[1];

  file.close();

  return {conf, seq};
}

void DNA_Analyzer::InitializeWeights() {
  for (int i = 0; i < 4; ++i) {
    m_weights[m_ABC[i]] = i;
  }
}

std::vector<std::vector<int> > DNA_Analyzer::InitScoreMatrix(const std::array<int, 3> &conf, const std::size_t rows, const std::size_t cols) {
  std::vector<std::vector<int> > mat(rows, std::vector<int>(cols, 0));
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) {
      mat[i][j] = ((i == j) ? conf[0] : conf[1]);
    }
  }
  return mat;
}

std::vector<std::vector<int> > DNA_Analyzer::InitAlighmentMatrix(const std::array<int, 3> &conf, const std::array<std::string, 2> &seq, const std::vector<std::vector<int> > &score_mat) {
  const int rows = seq[0].length() + 1;
  const int cols = seq[1].length() + 1;
  std::vector<std::vector<int> > mat(rows, std::vector<int>(cols, 0));

  for (int i = 0; i < rows; ++i) {
    mat[i][0] = i * conf[2];
  }

  for (int j = 0; j < cols; ++j) {
    mat[0][j] = j * conf[2];
  }

  for (int i = 1; i < rows; ++i) {
    const char c1 = seq[0][i - 1];
    for (int j = 1; j < cols; ++j) {
      const char c2 = seq[1][j - 1];

      const int score1 = mat[i - 1][j - 1] + score_mat[m_weights[c1]][m_weights[c2]];
      const int score2 = mat[i - 1][j] + conf[2];
      const int score3 = mat[i][j - 1] + conf[2];

      mat[i][j] = std::max({score1, score2, score3});
    }
  }

  return mat;
}

std::array<std::string, 3> DNA_Analyzer::ComputeAlignedSequences(const std::array<int, 3> &conf, const std::array<std::string, 2> &seq, const std::vector<std::vector<int> > &mat, const std::vector<std::vector<int> > &score_mat) {
  std::array<std::string, 3> alig;

  int i = seq[0].length(), j = seq[1].length();
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && mat[i][j] == mat[i - 1][j - 1] + score_mat[m_weights[seq[0][i - 1]]][m_weights[seq[1][j - 1]]]) {
      alig[0] += seq[0][--i];
      alig[2] += seq[1][--j];
    } else if (i > 0 && mat[i][j] == mat[i - 1][j] + conf[2]) {
      alig[0] += seq[0][--i];
      alig[2] += '-';
    } else {
      alig[0] += '-';
      alig[2] += seq[1][--j];
    }
    alig[1] += alig[0].back() == alig[2].back() ? '|' : ' ';
  }

  std::reverse(alig[0].begin(), alig[0].end());
  std::reverse(alig[2].begin(), alig[2].end());
  std::reverse(alig[1].begin(), alig[1].end());

  return alig;
}

void DNA_Analyzer::PrintNWResult(const std::array<std::string, 3> &alig, const int &score) {
  std::cout << "Max score: " << score << std::endl;

  std::cout << alig[0] << std::endl;
  std::cout << alig[1] << std::endl;
  std::cout << alig[2] << std::endl;
}

void DNA_Analyzer::RegexAlgorithm(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::string seq, pattern;
  file >> seq >> pattern;

  const std::size_t n = seq.size();
  const std::size_t m = pattern.size();

  std::vector<std::vector<bool>> dp(n + 1, std::vector<bool>(m + 1, false));
  dp[0][0] = true;

  for (auto i = seq.begin(); i != seq.end(); ++i) {
    const auto x = std::distance(seq.begin(), i) + 1;
    for (auto j = pattern.begin(); j != pattern.end(); ++j) {
      const auto y = std::distance(pattern.begin(), j) + 1;

      if (symbol(*i) == symbol(*j) || *j == '.') {
        dp[x][y] = dp[x - 1][y - 1];
      } else if (*j == '?') {
        dp[x][y] = dp[x - 1][y - 1] || dp[x][y - 1];
      } else if (*j == '*') {
        dp[x][y] = dp[x - 1][y - 1] || dp[x - 1][y];
      } else if (*j == '+') {
        dp[x][y] = (dp[x][y - 1] || dp[x - 1][y]) && (symbol(*i) == symbol(*(j - 1)) || *(j - 1) == '.');
      }
    }
  }

  std::cout << (dp[n][m] ? "True" : "False") << std::endl;
}

char DNA_Analyzer::symbol(char c) {
  return std::find(m_ABC.begin(), m_ABC.end(), c) != m_ABC.end() ? c : '\0';
}

void DNA_Analyzer::KSimilarAlgorithm(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::string s, p;
  file >> s >> p;
  std::cout << kSimilarity(s, p) << std::endl;
}

int DNA_Analyzer::kSimilarity(std::string s1, std::string s2) {
  if (s1.size() != s2.size()) return -1;
  std::unordered_set<std::string> visited;
  std::queue<std::pair<std::string, int>> q;
  q.push({s1, 0});
  visited.insert(s1);
  while (!q.empty()) {
    std::string curr = q.front().first;
    int swaps = q.front().second;
    q.pop();
    if (curr == s2) return swaps;
    int i = 0;
    while (curr[i] == s2[i]) i++;
    for (int j = i + 1; j < curr.size(); j++) {
      if (curr[j] == s2[i] && curr[j] != s2[j]) {
        std::swap(curr[i], curr[j]);
        if (visited.count(curr) == 0) {
          q.push({curr, swaps + 1});
          visited.insert(curr);
        }
        std::swap(curr[i], curr[j]);
      }
    }
  }
  return -1;
}

void DNA_Analyzer::WindowAlgorithm(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::string s, p;
  file >> s >> p;
  std::cout << minWindowSubstring(s, p) << std::endl;
}

std::string DNA_Analyzer::minWindowSubstring(std::string s, std::string t) {
  std::unordered_map<char, int> mp;
  for (char c : t) mp[c]++;
  int cnt = mp.size(), left = 0, right = 0, ansL = -1, ansR = -1,
      minLen = INT_MAX;
  while (right < s.size()) {
    if (mp.count(s[right])) {
      mp[s[right]]--;
      if (mp[s[right]] == 0) cnt--;
    }
    right++;

    while (cnt == 0) {
      if (right - left < minLen) {
        minLen = right - left;
        ansL = left;
        ansR = right;
      }
      if (mp.count(s[left])) {
        mp[s[left]]++;
        if (mp[s[left]] > 0) cnt++;
      }
      left++;
    }
  }
  if (ansL == -1)
    return "";
  else
    return s.substr(ansL, ansR - ansL);
}