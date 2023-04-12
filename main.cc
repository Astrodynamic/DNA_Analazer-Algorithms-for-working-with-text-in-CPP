#include <locale>

#include "interface.h"

int main() {
  std::locale::global(std::locale(""));
  Interface interface;
  interface.Exec();
  return 0;
}
