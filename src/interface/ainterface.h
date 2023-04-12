#pragma once

#include <cstdint>
#include <filesystem>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

class AbstractInterface {
 public:
  AbstractInterface() = default;
  explicit AbstractInterface(const AbstractInterface &other) = delete;
  explicit AbstractInterface(AbstractInterface &&other) = delete;
  AbstractInterface &operator=(const AbstractInterface &other) = delete;
  AbstractInterface &operator=(AbstractInterface &&other) = delete;
  virtual ~AbstractInterface() = default;

  virtual void Exec() = 0;

 protected:
  const bool RunMenu(const std::vector<std::function<bool(void)>> &func,
                     std::size_t menu);
  const std::size_t ShowMenu(const std::string &menu,
                             const std::size_t items = 0U);

  [[nodiscard]] const std::int64_t CheckInputItem(const std::int64_t min,
                                                  const std::int64_t max);
  [[nodiscard]] std::pair<bool, std::filesystem::path> CheckInputPathFile();
  [[nodiscard]] const bool Exit();

  static const std::vector<std::string> m_menus;
  std::vector<std::vector<std::function<bool(void)>>> m_funcs;
};
