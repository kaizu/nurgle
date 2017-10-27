#include <iostream>

#include <nurgle/csv.hpp>

int main(int argc, char* argv[])
{
    using namespace nurgle;

    std::string const pathto(argc > 0 ? argv[1] : ".");

    std::string const filename("/compounds.csv");
    std::vector<std::tuple<std::string, double>> const compounds = utils::csv<std::string, double>::read(pathto + filename);
    for (auto elem : compounds)
    {
        std::cout << std::get<0>(elem) << " : " << std::get<1>(elem) << std::endl;
    }

    return 0;
}
