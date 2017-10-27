#include <iostream>

#include <nurgle/csv.hpp>

int main(int argc, char* argv[])
{
    using namespace nurgle;

    std::string const sep("/");
    std::string const pathto(argc > 0 ? argv[1] : ".");

    std::string const filename1("compounds.csv");
    std::vector<std::tuple<std::string, double>> const compounds
        = utils::csv<std::string, double>::read(pathto + sep + filename1);
    for (auto elem : compounds)
    {
        std::cout << std::get<0>(elem) << " : " << std::get<1>(elem) << std::endl;
    }

    std::string const filename2("metabolism.csv");
    std::vector<std::tuple<std::string, std::string, double>> const metabolism
        = utils::csv<std::string, std::string, double>::read(pathto + sep + filename2);
    for (auto elem : metabolism)
    {
        std::cout << std::get<0>(elem) << " : " << std::get<1>(elem) << " : " << std::get<2>(elem) << std::endl;
    }

    return 0;
}
