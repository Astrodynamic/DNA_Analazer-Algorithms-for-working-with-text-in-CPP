#pragma once

#include "model.h"
#include "ainterface.h"

class Interface final : virtual public AbstractInterface {
 public:
  Interface();
  explicit Interface(const Interface &other) = delete;
  explicit Interface(Interface &&other) = delete;
  Interface &operator=(const Interface &other) = delete;
  Interface &operator=(Interface &&other) = delete;
  ~Interface();

  virtual void Exec() final override;

 private:
  enum MenuFuncs : std::size_t {
    kMainFuncMenu = 0U,
    kRabinKarpMenu,
    kNWAlgorithmMenu,
    kRegexMenu,
    kSimilarMenu,
    kWindowMenu,
    kMenuFuncsAll
  };

  enum MenuItem : std::size_t {
    kIntroduction = 0U,
    kMainMenu,
    kLoadMenu,
    kNotExistMenus,
    kCompletion
  };

  ModelAlgorithm m_analyzer;

  void InitFuncMenus();
  bool RunProcessFile1arg(std::function<void(const std::filesystem::path &)> func);
  bool RunProcessFile2arg(std::function<void(const std::filesystem::path &, const std::filesystem::path &)> func);
};