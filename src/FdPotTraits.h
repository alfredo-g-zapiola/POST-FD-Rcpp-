#include <functional>
#include <vector>

namespace fdpot_traits{
  using FDataType = std::vector<std::function<double(double)>>;

}
//struct FdPotTraits{
//  using FDataContainer = std::vector<std::function<double(double)>>;
//};