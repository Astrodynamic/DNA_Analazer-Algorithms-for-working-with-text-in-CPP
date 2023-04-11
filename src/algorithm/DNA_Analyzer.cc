#include "DNA_Analyzer.h"

void DNA_Analyzer::RabinKarpAlgorithm(const std::filesystem::path& path_1, const std::filesystem::path& path_2) {
  std::string text = std::move(ReadFile(path_1));
  std::string samp = std::move(ReadFile(path_2));

  const std::size_t text_size = text.size();
  const std::size_t samp_size = samp.size();

  std::vector<std::size_t> exp(std::max(text_size, samp_size), 1);
  std::transform(exp.begin(), exp.end() - 1, exp.begin() + 1, [this](std::size_t &item) { return (item * rk_base) % rk_mod; });

  std::vector<std::size_t> hash_text(text_size + 1, 0);
  CalculateMassHash(text, hash_text, text_size, exp);

  std::vector<std::size_t> hash_samp(samp_size + 1, 0);
  CalculateMassHash(samp, hash_samp, samp_size, exp);

  for (std::size_t idx = 0; idx + samp_size - 1 < text_size; ++idx) {
    std::size_t hash_curr = (hash_text[idx + samp_size] + rk_mod - hash_text[idx]) % rk_mod;
    if ((hash_curr == (hash_samp[samp_size] * exp[idx] % rk_mod)) && text.substr(idx, samp_size) == samp) {
      std::cout << idx << " ";
    }
  }
  std::cout << std::endl;
}

void DNA_Analyzer::CalculateMassHash(const std::string &src, std::vector<std::size_t> &hash, const std::size_t size, const std::vector<std::size_t> &exp) {
  for (std::size_t i = 0; i < size; ++i) {
    hash[i + 1] = (hash[i] + src[i] * exp[i]) % rk_mod;
  }
}

std::string DNA_Analyzer::ReadFile(const std::filesystem::path& path) {
  std::ifstream file(path);
  return std::string((std::istreambuf_iterator<char>(file)), {});
}

void DNA_Analyzer::NWAlgorithm(const std::filesystem::path& path) {
  std::ifstream file(path);

  Score score(0, 0, 0);
  file >> score.match >> score.mismatch >> score.gap;
  std::string a, b;
  file >> a >> b;

  std::cout << optimalScore(score, a, b);

  auto k = optimalAlignment(score, a, b);
  std::cout << k.first << std::endl;
  std::cout << k.second.first << std::endl;
  std::cout << k.second.second << std::endl;
}

int DNA_Analyzer::optimalScore(Score score, std::string s1, std::string s2) {
  int n = s1.size(), m = s2.size();
  std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1));
  for (int i = 1; i <= n; i++) dp[i][0] = dp[i - 1][0] + score.gap;
  for (int i = 1; i <= m; i++) dp[0][i] = dp[0][i - 1] + score.gap;
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      int match = dp[i - 1][j - 1] +
                  (s1[i - 1] == s2[j - 1] ? score.match : score.mismatch);
      int ins = dp[i][j - 1] + score.gap;
      int del = dp[i - 1][j] + score.gap;
      dp[i][j] = std::max({match, ins, del});
    }
  }
  return dp[n][m];
}

std::pair<int, std::pair<std::string, std::string>>
DNA_Analyzer::optimalAlignment(Score score, std::string s1, std::string s2) {
  int n = s1.size(), m = s2.size();
  std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1));
  std::vector<std::vector<int>> p(
      n + 1, std::vector<int>(m + 1));  // 0-diag, 1-up, 2-left

  for (int i = 1; i <= n; i++) dp[i][0] = dp[i - 1][0] + score.gap, p[i][0] = 1;
  for (int i = 1; i <= m; i++) dp[0][i] = dp[0][i - 1] + score.gap, p[0][i] = 2;
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      int match = dp[i - 1][j - 1] +
                  (s1[i - 1] == s2[j - 1] ? score.match : score.mismatch);
      int ins = dp[i][j - 1] + score.gap;
      int del = dp[i - 1][j] + score.gap;
      if (match >= ins && match >= del) {
        dp[i][j] = match;
        p[i][j] = 0;
      } else if (ins >= match && ins >= del) {
        dp[i][j] = ins;
        p[i][j] = 2;
      } else {
        dp[i][j] = del;
        p[i][j] = 1;
      }
    }
  }

  int i = n, j = m;
  std::string res1, res2;
  while (i > 0 || j > 0) {
    if (p[i][j] == 0) {
      res1 += s1[i - 1];
      res2 += s2[j - 1];
      i--;
      j--;
    } else if (p[i][j] == 1) {
      res1 += s1[i - 1];
      res2 += '-';
      i--;
    } else {
      res1 += '-';
      res2 += s2[j - 1];
      j--;
    }
  }
  reverse(res1.begin(), res1.end());
  reverse(res2.begin(), res2.end());

  return {dp[n][m], {res1, res2}};
}

void DNA_Analyzer::RegexAlgorithm(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::string s, p;
  file >> s >> p;
  std::cout << isMatch(s, p) << std::endl;
}

bool DNA_Analyzer::isMatch(std::string s, std::string p) {
  int m = p.size(), n = s.size();
  std::vector<std::vector<bool>> dp(m + 1, std::vector<bool>(n + 1));
  dp[0][0] = true;
  for (int i = 1; i <= m; i++) {
    if (p[i - 1] == '*') dp[i][0] = dp[i - 2][0];
  }
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= n; j++) {
      if (s[j - 1] == p[i - 1] || p[i - 1] == '.')
        dp[i][j] = dp[i - 1][j - 1];
      else if (p[i - 1] == '*') {
        dp[i][j] = dp[i - 2][j] || (dp[i][j - 1] && (s[j - 1] == p[i - 2] || p[i - 2] == '.'));
      } else if (p[i - 1] == '?') {
        dp[i][j] = dp[i - 1][j - 1] || dp[i - 1][j] || dp[i][j - 1];
      }
    }
  }
  return dp[m][n];
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