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

 private:
  std::vector<int> RabinKarpSearch(std::string a, std::string b) {
    int n = a.size(), m = b.size();
    std::vector<int> res;
    if (m > n) return res;

    long long h = 0, p = 1;
    const int base = 101;
    for (int i = 0; i < m; i++) {
      h = h * base + b[i];
      p *= base;
    }

    long long cur = 0;
    for (int i = 0; i < n; i++) {
      cur = cur * base + a[i];
      if (i >= m) cur -= a[i - m] * p;
      if (i >= m - 1 && cur == h) res.push_back(i - m + 1);
    }
    return res;
  }
};

// // Part 2: NW sequence alignment project

// struct Score {
//   int match, mismatch, gap;
//   Score(int _match, int _mismatch, int _gap)
//       : match(_match), mismatch(_mismatch), gap(_gap) {}
// };

// class NeedlemanWunsch {
//  public:
//   int optimalScore(Score score, string s1, string s2) {
//     int n = s1.size(), m = s2.size();
//     vector<vector<int>> dp(n + 1, vector<int>(m + 1));
//     for (int i = 1; i <= n; i++) dp[i][0] = dp[i - 1][0] + score.gap;
//     for (int i = 1; i <= m; i++) dp[0][i] = dp[0][i - 1] + score.gap;
//     for (int i = 1; i <= n; i++) {
//       for (int j = 1; j <= m; j++) {
//         int match = dp[i - 1][j - 1] +
//                     (s1[i - 1] == s2[j - 1] ? score.match : score.mismatch);
//         int ins = dp[i][j - 1] + score.gap;
//         int del = dp[i - 1][j] + score.gap;
//         dp[i][j] = max({match, ins, del});
//       }
//     }
//     return dp[n][m];
//   }

//   pair<int, pair<string, string>> optimalAlignment(Score score, string s1,
//                                                    string s2) {
//     int n = s1.size(), m = s2.size();
//     vector<vector<int>> dp(n + 1, vector<int>(m + 1));
//     vector<vector<int>> p(n + 1, vector<int>(m + 1));  // 0-diag, 1-up,
//     2-left

//     for (int i = 1; i <= n; i++)
//       dp[i][0] = dp[i - 1][0] + score.gap, p[i][0] = 1;
//     for (int i = 1; i <= m; i++)
//       dp[0][i] = dp[0][i - 1] + score.gap, p[0][i] = 2;
//     for (int i = 1; i <= n; i++) {
//       for (int j = 1; j <= m; j++) {
//         int match = dp[i - 1][j - 1] +
//                     (s1[i - 1] == s2[j - 1] ? score.match : score.mismatch);
//         int ins = dp[i][j - 1] + score.gap;
//         int del = dp[i - 1][j] + score.gap;
//         if (match >= ins && match >= del) {
//           dp[i][j] = match;
//           p[i][j] = 0;
//         } else if (ins >= match && ins >= del) {
//           dp[i][j] = ins;
//           p[i][j] = 2;
//         } else {
//           dp[i][j] = del;
//           p[i][j] = 1;
//         }
//       }
//     }

//     int i = n, j = m;
//     string res1, res2;
//     while (i > 0 || j > 0) {
//       if (p[i][j] == 0) {
//         res1 += s1[i - 1];
//         res2 += s2[j - 1];
//         i--;
//         j--;
//       } else if (p[i][j] == 1) {
//         res1 += s1[i - 1];
//         res2 += '-';
//         i--;
//       } else {
//         res1 += '-';
//         res2 += s2[j - 1];
//         j--;
//       }
//     }
//     reverse(res1.begin(), res1.end());
//     reverse(res2.begin(), res2.end());

//     return {dp[n][m], {res1, res2}};
//   }
// };

// // Part 3: Matching regular expressions
// class RegexMatching {
//  public:
//   bool isMatch(string s, string p) {
//     int m = p.size(), n = s.size();
//     vector<vector<bool>> dp(m + 1, vector<bool>(n + 1));
//     dp[0][0] = true;
//     for (int i = 1; i <= m; i++) {
//       if (p[i - 1] == '*') dp[i][0] = dp[i - 2][0];
//     }
//     for (int i = 1; i <= m; i++) {
//       for (int j = 1; j <= n; j++) {
//         if (s[j - 1] == p[i - 1] || p[i - 1] == '.')
//           dp[i][j] = dp[i - 1][j - 1];
//         else if (p[i - 1] == '*') {
//           dp[i][j] =
//               dp[i - 2][j] ||
//               (dp[i][j - 1] && (s[j - 1] == p[i - 2] || p[i - 2] == '.'));
//         } else if (p[i - 1] == '?') {
//           dp[i][j] = dp[i - 1][j - 1] || dp[i - 1][j] || dp[i][j - 1];
//         }
//       }
//     }
//     return dp[m][n];
//   }
// };

// // Part 4: K-similar strings
// class KSimilarity {
//  public:
//   int kSimilarity(string s1, string s2) {
//     if (s1.size() != s2.size()) return -1;
//     unordered_set<string> visited;
//     queue<pair<string, int>> q;
//     q.push({s1, 0});
//     visited.insert(s1);
//     while (!q.empty()) {
//       string curr = q.front().first;
//       int swaps = q.front().second;
//       q.pop();
//       if (curr == s2) return swaps;
//       int i = 0;
//       while (curr[i] == s2[i]) i++;
//       for (int j = i + 1; j < curr.size(); j++) {
//         if (curr[j] == s2[i] && curr[j] != s2[j]) {
//           swap(curr[i], curr[j]);
//           if (visited.count(curr) == 0) {
//             q.push({curr, swaps + 1});
//             visited.insert(curr);
//           }
//           swap(curr[i], curr[j]);
//         }
//       }
//     }
//     return -1;
//   }
// };

// // Part 5: Minimum window substring
// class MinWindowSubstring {
//  public:
//   string minWindowSubstring(string s, string t) {
//     unordered_map<char, int> mp;
//     for (char c : t) mp[c]++;
//     int cnt = mp.size(), left = 0, right = 0, ansL = -1, ansR = -1,
//         minLen = INT_MAX;
//     while (right < s.size()) {
//       if (mp.count(s[right])) {
//         mp[s[right]]--;
//         if (mp[s[right]] == 0) cnt--;
//       }
//       right++;

//       while (cnt == 0) {
//         if (right - left < minLen) {
//           minLen = right - left;
//           ansL = left;
//           ansR = right;
//         }
//         if (mp.count(s[left])) {
//           mp[s[left]]++;
//           if (mp[s[left]] > 0) cnt++;
//         }
//         left++;
//       }
//     }
//     if (ansL == -1)
//       return "";
//     else
//       return s.substr(ansL, ansR - ansL);
//   }
// };
