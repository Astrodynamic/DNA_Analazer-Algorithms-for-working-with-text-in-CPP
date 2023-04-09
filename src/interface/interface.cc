#include "interface.h"

const std::vector<std::string> Interface::AbstractInterface::m_menus{
    " -------------------------------------------------------------- \n"
    "|                       DNA Analyzer 1.0                       |\n"
    " -------------------------------------------------------------- \n",
    " -------------------------------------------------------------- \n"
    "|                       Select menu item                       |\n"
    " -------------------------------------------------------------- \n"
    "| 0. Exit                                                      |\n"
    "| 1. The Exact DNA search                                      |\n"
    "| 2. The NW sequence alignment                                 |\n"
    "| 3. Matching regular expressions                              |\n"
    "| 4. K-similar strings                                         |\n"
    "| 5. Minimum window substring                                  |\n"
    " -------------------------------------------------------------- \n"
    " > ",
    " -------------------------------------------------------------- \n"
    "|                       Select menu item                       |\n"
    " -------------------------------------------------------------- \n"
    "| 0. Exit                                                      |\n"
    "| 1. Enter the full path to the file ...                       |\n"
    " -------------------------------------------------------------- \n"
    " > ",
    " -------------------------------------------------------------- \n"
    "|             A file with that name does not exist             |\n"
    " -------------------------------------------------------------- \n",
    " -------------------------------------------------------------- \n"
    "|            Successful completion of the programme            |\n"
    " -------------------------------------------------------------- \n"};

Interface::Interface() {
  InitFuncMenus();
  ShowMenu(m_menus[MenuItem::kIntroduction]);
}

void Interface::InitFuncMenus() {
  m_funcs.resize(MenuFuncs::kMenuFuncsAll);

  m_funcs[MenuFuncs::kMainFuncMenu] = {
    std::bind(&Interface::Exit, this),
    std::bind(&Interface::RunMenu, this, std::ref(m_funcs[MenuFuncs::kRabinKarpMenu]), MenuItem::kLoadMenu)
  };

  m_funcs[MenuFuncs::kRabinKarpMenu] = {
    std::bind(&Interface::Exit, this),
    [this]() -> const bool {
      return RunProcessFile(std::bind(&DNA_Analyzer::RabinKarpAlgorithm, std::ref(m_analyzer), std::placeholders::_1, std::placeholders::_2));
    }
  };
}

const bool Interface::RunProcessFile(std::function<void(const std::filesystem::path& path)> func) {
  const auto [success, path] {CheckInputPathFile()};
  if (success) {
    func(path);
  } else {
    ShowMenu(m_menus[MenuItem::kNotExistMenus]);
  }
  return !success;
}

const bool Interface::RunProcessFile(std::function<void(const std::filesystem::path&, const std::filesystem::path&)> func) {
  const auto [success_1, path_1] {CheckInputPathFile()};
  const auto [success_2, path_2] {CheckInputPathFile()};
  if (success_1 && success_2) {
    func(path_1, path_2);
  } else {
    ShowMenu(m_menus[MenuItem::kNotExistMenus]);
  }
  return !success_1 || !success_2;
}

Interface::~Interface() {
  ShowMenu(m_menus[MenuItem::kCompletion]);
}

void Interface::Exec() {
  RunMenu(m_funcs[MenuFuncs::kMainFuncMenu], MenuItem::kMainMenu);
}