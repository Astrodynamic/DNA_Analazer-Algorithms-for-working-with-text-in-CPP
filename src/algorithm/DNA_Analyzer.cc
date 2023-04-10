#include "DNA_Analyzer.h"

void DNA_Analyzer::RabinKarpAlgorithm(const std::filesystem::path& path_1, const std::filesystem::path& path_2) {
  std::ifstream file_a(path_1), file_b(path_2);

  std::string a((std::istreambuf_iterator<char>(file_a)), std::istreambuf_iterator<char>());
  std::string b((std::istreambuf_iterator<char>(file_b)), std::istreambuf_iterator<char>());

  const int p = 31; // простое число
  const int m = 1e9 + 9; // модуль
  int S = a.size(), T = b.size();
  
  std::vector<long long> p_pow(std::max(S, T)); // предподсчет степеней числа p
  p_pow[0] = 1;
  for (int i = 1; i < (int)p_pow.size(); i++) {
    p_pow[i] = (p_pow[i - 1] * p) % m;
  }
  
  std::vector<long long> h(S + 1, 0); // хеши от всех префиксов строки text
  for (int i = 0; i < S; i++) {
    int code;
    switch(a[i]) { // перевод из алфавита {A, C, G, T} в числа
        case 'A': code = 1; break;
        case 'C': code = 2; break;
        case 'G': code = 3; break;
        case 'T': code = 4; break;
    }
    h[i + 1] = (h[i] + code * p_pow[i]) % m;
  }
  
  long long h_s = 0; // хеш подстроки длины T строки text
  for (int i = 0; i < T; i++) {
    int code;
    switch(b[i]) { // перевод из алфавита {A, C, G, T} в числа
        case 'A': code = 1; break;
        case 'C': code = 2; break;
        case 'G': code = 3; break;
        case 'T': code = 4; break;
    }
    h_s = (h_s + code * p_pow[i]) % m;
  }
  
  std::vector<int> positions; // позиции вхождения подстроки pattern в строку text
  for (int i = 0; i + T - 1 < S; i++) {
    long long cur_h = (h[i + T] + m - h[i]) % m; // вычисление хеша для подстроки
    if (cur_h == h_s * p_pow[i] % m) {
      positions.push_back(i);
    }
  }

  for (auto &it: positions) {
    std::cout << it << " ";
  }
  std::cout << std::endl;
}
