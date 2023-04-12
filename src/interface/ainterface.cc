#include "ainterface.h"

const bool AbstractInterface::RunMenu(
    const std::vector<std::function<bool(void)>> &func, std::size_t menu) {
  bool flag{};
  std::size_t item{};
  std::size_t items{func.size()};
  do {
    item = ShowMenu(m_menus[menu], items);
  } while (items && (flag = func[item]()));

  return !flag;
}

const std::size_t AbstractInterface::ShowMenu(const std::string &menu,
                                              const std::size_t items) {
  std::cout << menu;
  return items ? CheckInputItem(-1, items) : 0;
}

[[nodiscard]] const std::int64_t AbstractInterface::CheckInputItem(
    const std::int64_t min, const std::int64_t max) {
  std::string line;
  std::getline(std::cin, line);

  std::int64_t result;
  while (!sscanf(line.c_str(), "%lld", &result) || result <= min ||
         result >= max) {
    std::cout << "Incorrect input, try again: ";
    std::getline(std::cin, line);
  }

  return result;
}

[[nodiscard]] std::pair<bool, std::filesystem::path>
AbstractInterface::CheckInputPathFile() {
  std::string line;
  std::getline(std::cin, line);
  std::filesystem::path path(line);
  return std::make_pair(std::filesystem::exists(path), path);
}

[[nodiscard]] const bool AbstractInterface::Exit() { return false; }